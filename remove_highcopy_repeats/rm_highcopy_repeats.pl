use strict;
use warnings;
use 5.010;

@ARGV == 3 or &print_usage();

my ($sets_file, $wgac_file, $out_file) = @ARGV;
open IN,$sets_file or die "Cannot open file $sets_file: $!";
$_=<IN>;
say "Reading $sets_file";
my @rm;
my $num_highcopy=0;
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	if($c[0] >= 3 or $c[2] >= 50)
	{
		$num_highcopy ++;
		my %h;
		for(@c[3..$#c])
		{
			my @chr_loc = split / /;
			$chr_loc[0] =~ s/://;
			push @{$h{$chr_loc[0]}},@chr_loc[1..$#chr_loc];
		}
		
		push @rm,{%h};
	}
}
close IN;
say "# of high-copy repeats: $num_highcopy";

say "Removing high-copy repeats in $wgac_file";
my $num_removed_hits = 0;
open IN,$wgac_file or die "Cannot open file $wgac_file: $!";
open OUT,'>'.$out_file;
my $line;
my $ln_nmbr = 0;
while($line = <IN>)
{
	$ln_nmbr ++;
	$ln_nmbr % 10000 == 0 and say "Scanned $ln_nmbr lines";
	my @c = split /\t/,$line;
	
	$c[0] =~ s/chrom//;
	$c[0] =~ s/\.fa//;
	
	$c[3] =~ s/chrom//;
	$c[3] =~ s/\.fa//;
	
	my $f=0;
	for my $sd (@rm)
	{
		my %h = %$sd;
		
		if(defined $h{$c[0]} and defined $h{$c[3]})
		{
			my $f1 = &chk($c[1],$c[2],@{$h{$c[0]}});
			$c[4] > $c[5] and @c[4,5]=@c[5,4];
			my $f2 = &chk($c[4],$c[5],@{$h{$c[3]}});
			
			if($f1 and $f2)
			{
				$f = 1;
				last;
			}
		}
		else
		{
			next;
		}
	}
	
	if($f)
	{
		$num_removed_hits ++;
	}
	else
	{
		print OUT $line;
	}
}
say "Number of lines removed in $wgac_file: $num_removed_hits";

# Compare locations of a given hit to high-copy repeats 
sub chk
{
	my $ps = shift;
	my $pe = shift;

	for(@_)
	{
		my ($s,$e) = split /\-/;
		if(abs($ps-$s) < 1000 or abs($pe-$e) < 1000)
		{
			return 1;
		}
	}
	return 0;
}

sub print_usage()
{
	print "Run the program: \n";
	print "  perl rm_highcopy_repeats.pl argv1 argv2 argv3 \n";
	print "  The first argument should be the output file of clean_sets_of_similar_hits.pl \n";
	print "  The second argument should be the FINAL output file of the WGAC pipeline \n";
	print "  The third argument should be the output file you want to specify. \n";
	exit;
}
