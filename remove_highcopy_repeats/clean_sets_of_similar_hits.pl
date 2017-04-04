# clean_sets_of_similar_hits.pl
# This program is used to clean sets of similar hits from get_sets_of_similar_hits.pl
# 

use strict;
use warnings;
use 5.010;

@ARGV == 2 or &print_usage();

my @ss;
open IN,$ARGV[0] or die "Cannot open file $ARGV[0]: $!";    # the output file of get_sets_of_similar_hits.pl
$_=<IN>;
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	my %h;
	$h{'numChrom'} = $c[0];
	$h{'numPair'} = $c[1]/2;   # Note this is not really the number of pairs of hits.
	for (@c[2..$#c])
	{
		my @chr_loc = split / /;
		$chr_loc[0] =~ s/\:$//;
		$h{$chr_loc[0]} = [&rm_repeat(@chr_loc[1..$#chr_loc])];
	}
	
	push @ss, {%h};
}
say "Finished removing duplicate hits.";


ILOOP: for(my $i=0; $i < @ss; $i++)
{
	my %hi=%{$ss[$i]};
	JLOOP: for(my $j=$i+1; $j < @ss; $j++)
	{
		my %hj=%{$ss[$j]};
		for my $ikey (keys %hi)
		{
			next if ($ikey =~ /^num/);
			if(defined $hj{$ikey})
			{
				if(&cmpr(@{$hj{$ikey}}, @{$hi{$ikey}}) == 1)
				{
					$ss[$j]->{'numPair'} += $ss[$i]->{'numPair'};  # This stat will not used.
					for(keys %hi)
					{
						next if (/^num/);
						if(defined $hj{$_})
						{
							$ss[$j]->{$_} = [&join_loc(@{$hj{$_}}, @{$hi{$_}})];
						}
						else
						{
							$ss[$j]->{$_} = $hi{$_};
							$ss[$j]->{'numChrom'} ++;
						}
					}
					say "SET $i has been merged to SET $j.";
					splice (@ss,$i,1);
					$i--;
					next ILOOP;
				}
			}
		}
	}
}
say "Finished merging sets of similar hits.";
open OUT,'>'.$ARGV[1];
say OUT "NumChrom\tNumPairs\tNumCopies\tHits(Chrom: Pos)";
for my $sd_ref (@ss)
{
	my %h = %$sd_ref;
	
	my $sd_cnt;
	my $known_chrom_cnt;
	for (keys %h)
	{
		next if (/^num/);
		if(/^[\dXY]/) {$known_chrom_cnt ++;}      # only count known chromosomes, i.e. autosomes and allosomes
		$sd_cnt += @{$h{$_}};
	}

	print OUT $known_chrom_cnt,"\t",$h{'numPair'},"\t",$sd_cnt,"\t";
	
	for my $chrom (sort keys %h)
	{
		next if ($chrom =~ /^num/);
		print OUT "$chrom:";
		for (@{$sd_ref->{$chrom}})
		{
			print OUT " ";
			print OUT join "-",@$_;
		}
		print OUT "\t";
	}
	print OUT "\n";
}
say "Finished printing repeats into output file.";

# remove a hit if it has very similar positions (+-1k) to a previous one.
sub rm_repeat {
	my @in = @_;
	my @out=();
	for(@in)
	{
		my ($s,$e)=split /\-/;
		if(@out)
		{
			if(abs($out[-1][0]-$s) < 1000 and abs($out[-1][1]-$e) <1000)
			{
				next;
			}
			else
			{
				push @out,[$s,$e];
			}
		}
		else
		{
			push @out,[$s,$e];
		}
	}
	return @out;
}

# compare if two sets of similar hits are of the same one
sub cmpr {
	my @in = sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @_;
	my $out = 0;
	
	for(my $i=0; $i<@in-1; $i++)
	{
		if(abs($in[$i][0]-$in[$i+1][0]) < 1000 and abs($in[$i][1]-$in[$i+1][1]) < 1000)
		{
			$out = 1;
			last;
		}
	}
	return $out;
}

# join hits
sub join_loc {
	my @in = sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @_;
	my @out=();
	for(@in)
	{
		my ($s,$e)=@$_;
		if(@out)
		{
			if(abs($out[-1][0]-$s) < 1000 and abs($out[-1][1]-$e) <1000)
			{
				next;
			}
			else
			{
				push @out,[$s,$e];
			}
		}
		else
		{
			push @out,[$s,$e];
		}
	}
	return @out;
}

sub print_usage()
{
	print "Run the program: \n";
	print "  perl clean_sets_of_similar_hits.pl argv1 argv2 \n";
	print "  The first argument should be the output file of get_sets_of_similar_hits.pl \n";
	print "  The second argument should be the output file you want to specify. \n";
	exit;
}
