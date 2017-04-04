use warnings;
use strict;

our $p = "wgac_call2region.pl";

@ARGV == 3 or die "Run: $p *.wgac.clean lt|gt|eq identity
e.g., for joining wgac calls >0.98 identity into regions, the command should be like the following:
$p ss.wgac.clean gt 0.98\n";
open IN,$ARGV[0] or die "cannot find the input file $ARGV[0]: $!\n";
my $than_sign = $ARGV[1];
$than_sign eq 'lt' or $than_sign eq 'gt' or $than_sign eq 'eq' or die "The 2nd parameter must be lt for <, gt for > or eq for =.\n";
my $identity = $ARGV[2];
$identity =~ /^\d*\.\d+$/ or die "The 3rd parameter must be a positive float-point number";

my @wgac;
my $num_pair;
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	
	if($than_sign eq 'lt')
	{
		next if $c[-1] >= $identity;
	}
	elsif($than_sign eq 'gt')
	{
		next if $c[-1] <= $identity;
	}
	elsif($than_sign eq 'eq')
	{
		next if $c[-1] != $identity;
	}

	

	$num_pair ++;
	$c[0] =~ s/.fa$//;
	$c[0] =~ s/^chrom//;
	
	$c[3] =~ s/.fa$//;
	$c[3] =~ s/^chrom//;
	
	$c[4] > $c[5] and @c[4,5]=@c[5,4];
	push @wgac,([@c[0..2]],[@c[3..5]]);
}
print "According to filtering criteria, the number of pair alignments is: $num_pair\n";

if($num_pair == 0)
{
	print "Due to no pairwise alignment, no file is generated.\n";
	exit;
}

@wgac = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @wgac;

for(my $i=1; $i <@wgac; $i++)
{
	if($wgac[$i][0] eq $wgac[$i-1][0])
	{
		if($wgac[$i][1] <= $wgac[$i-1][2])
		{
			if($wgac[$i][2] > $wgac[$i-1][2])
			{
				$wgac[$i-1][2] = $wgac[$i][2];
			}
			splice (@wgac,$i,1);
			$i--;
		}
	}
}

print "Number of non-overlapping regions is: ", scalar(@wgac),"\n";

my $tl;
open OUT,'>'.$ARGV[0].".region.$than_sign-$identity";
for(@wgac)
{
	print OUT join "\t",@$_;
	$tl += $_->[2]-$_->[1]+1;
	print OUT "\t",$_->[2]-$_->[1]+1,"\n";
}

print "Length of non-overlapping regions is: ", $tl,"\n";

