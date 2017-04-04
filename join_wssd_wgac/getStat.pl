use warnings;
use strict;


@ARGV == 3 or die "Run: getStat.pl  WSSD  WGAC  Output\n    All ARGVS are files.\n";
open IN,$ARGV[0] or die "cannot find the input file $ARGV[0]: $!\n";
$_=<IN>;
my @wgac;
my $num_pair;
while(<IN>)
{
	chomp;
	my @c = split /\t/;

	$num_pair ++;
	$c[0] =~ s/^chr//;
	$c[1] ++;         # change to 1-origin coordinate as same as WGAC

	push @wgac,[@c[0..2]];
}

open IN,$ARGV[1] or die "cannot find the input file $ARGV[1]: $!\n";
while(<IN>)
{
	chomp;
	my @c = split /\t/;

	$num_pair ++;
	$c[0] =~ s/^chr//;

	push @wgac,[@c[0..2]];
}
print "Number of interval before joining is: $num_pair\n";
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

print "Number of non-overlapping regions after joining is: ", scalar(@wgac),"\n";

my $tl;
open OUT,'>'.$ARGV[2];
for(@wgac)
{
	print OUT join "\t",@$_;
	$tl += $_->[2]-$_->[1]+1;
	print OUT "\t",$_->[2]-$_->[1]+1,"\n";
}

print "Length of non-overlapping regions after joining is: ", $tl,"\n";

