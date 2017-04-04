use strict;
use warnings;

@ARGV == 2 or die "Run the program:\n  perl filter.pl  *wgac.clean  identity \n  ARG1 is a file and ARG2 is a float-point number.\n";
my ($wgac, $identity) = @ARGV;
$identity =~ /^\d*\.\d+$/ or die "The 2nd parameter must be a positive float-point number";

open IN, $wgac or die "cannot find the input file $wgac: $!\n";
my $outfile = "$wgac.ge-$identity";
open OUT,">$outfile";

while(<IN>)
{
	chomp;
	my @c = split /\t/;
	if($c[-1] >= 0.98)
	{
		print OUT $_,"\n";
	}
}
