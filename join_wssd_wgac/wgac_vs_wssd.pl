use strict;
use warnings;
use 5.010;

@ARGV == 3 or die "Run the program: perl wgac_vs_wssd.pl  WGAC  WSSD  output\n";

my %wgac;
my $num;
my $totlen;
open IN,$ARGV[0] or die "cannot find the input file $ARGV[0]: $!\n";
my $ind = 1;
while(<IN>)
{
	chomp;
	my @c=split /\t/;

	$c[0]=~s/^chrom//;
	$c[0]=~s/\.fa$//;
	$c[3]=~s/^chrom//;
	$c[3]=~s/\.fa$//;
	
	push @{$wgac{$c[0]}},[@c[1,2],$ind];
	$c[4] > $c[5] and @c[4,5]=@c[5,4];
	push @{$wgac{$c[3]}},[@c[4,5],$ind];
	$totlen += $c[2]-$c[1]+1;
	$totlen += $c[5]-$c[4]+1;
	$num += 2;
	$ind ++;
}
print "High-identity WGAC alignment intervals:\n\tNumber is $num, Length is $totlen\n\n";

for(keys %wgac)
{
	@{$wgac{$_}} = sort {$a->[0] <=> $b->[0]} @{$wgac{$_}};
}

my %wssd;
$num=0;
$totlen=0;
open IN,$ARGV[1] or die "cannot find the input file $ARGV[1]: $!\n";
$_=<IN>;
while(<IN>)
{
	chomp;
	my @c=split /\t/;
	next if ($c[3] eq 'del');
	next if ($c[2]-$c[1]<10000);
	$c[0] =~ s/chr//;
	push @{$wssd{$c[0]}},[$c[1]+1,$c[2]];       # change to 1-origin coordinate as same as WGAC
	$totlen += $c[2]-$c[1];
	$num++;
}
print "WSSD SD:\n\tNumber is $num, Length is $totlen\n\n";

for(keys %wssd)
{
	@{$wssd{$_}} = sort {$a->[0] <=> $b->[0]} @{$wssd{$_}};
}

say "Finding the WGAC hits overlapping WSSD hits.";
my @ss;
# col0: chr, col1: start & end of SD, col2: Overlapped Part of SD, col3: flag (WGAC/WSSD)
for my $c (keys %wgac)
{
	GLOOP: for(@{$wgac{$c}})                    # looping for wgac 
	{
		my ($s1,$e1,$idx1)=@$_;
		last GLOOP if (! defined $wssd{$c});
		KLOOP: for(my $i=0;$i<@{$wssd{$c}};$i++)  # looping for wssd
		{
			my ($s2,$e2)=@{$wssd{$c}->[$i]};
			if($e1 <= $s2)
			{
				next GLOOP;
			}
			elsif($s1 < $s2)
			{
				if($e1 < $e2)
				{
					push @ss,[$c,[$s1,$e1],[$s2,$e1],'WGAC',$idx1];
					next GLOOP;
				}
				else
				{
					push @ss,[$c,[$s1,$e1],[$s2,$e2],'WGAC',$idx1];
					next KLOOP;
				}
			}
			elsif($s1 < $e2)
			{
				if($e1 < $e2)
				{
					push @ss,[$c,[$s1,$e1],[$s1,$e1],'WGAC',$idx1];
					next GLOOP;
				}
				else
				{
					push @ss,[$c,[$s1,$e1],[$s1,$e2],'WGAC',$idx1];
					next KLOOP;
				}
			}
			else
			{
				splice @{$wssd{$c}},$i,1;
				$i--;
				next KLOOP;
			}
		}
	}
}

say "Calculating the proportion of WSSD-overlapping region in each WGAC hit";
my %st;
for(@ss)
{
	my ($c,$ref_sd,$ref_olap,$f,$idx)=@$_;
	$st{join(" ",($c,@$ref_sd,$idx))} += $ref_olap->[1] - $ref_olap->[0] + 1;
}
undef(@ss);

my $lmt = 0.1;
my %eff_ind;
for(keys %st)
{
	my $olap_len = $st{$_};
	my ($c,$s,$e,$idx)=split / /;
	my $l = $e-$s+1;
	if($l <= 2000)
	{
		$lmt = 0.2;
	}
	elsif($l <= 5000)
	{
		$lmt = 0.3;
	}
	elsif($l <= 10000)
	{
		$lmt = 0.4;
	}
	else
	{
		$lmt = 0.5;
	}
	next if $olap_len/$l < $lmt;
	$eff_ind{$idx} ++;
}
undef(%st);

say "Writing into output file: $ARGV[2]";
open OUT,'>'.$ARGV[2];
open IN,$ARGV[0] or die "cannot find the input file $ARGV[0]: $!\n";
$ind = 1;
while(<IN>)
{
	chomp;
	if(defined $eff_ind{$ind})           # and $eff_ind{$ind} == 2
	{
		print OUT $_,"\t",$eff_ind{$ind},"\n";
	}
	$ind ++;
}
