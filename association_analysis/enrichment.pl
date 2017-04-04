# Test if specific genomic sequences (like CNVs) are enriched in genomic features 
# (like gene regions, SDS) through simulation.

# Please contact Jicai Jiang (jicai.jiang@gmail.com) if you have any questions or concerns.

# all postions are 0-based.
# start and end are included in invervals.

use strict;
use warnings;

@ARGV == 5 or &print_usage();
my ($gf_file, $contig_file, $tested_file, $simu_num, $simu_out_prefix) = @ARGV;
$simu_num =~ /^\d+$/ or die "The 4th argv must be integer.\n";

# step 1
# read genome position of genomic features, e.g. SD, extron, intron

open IN,$gf_file or die "Cannot open file $gf_file: $!\n";
my %gf;
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	push @{$gf{$c[0]}},[@c[1,2]];
}
close IN;

for (keys %gf)
{
	@{$gf{$_}} = sort {$a->[0] <=> $b->[0]} @{$gf{$_}};
}

# step 2
# read genome excluding gap, i.e. contig position

my @contig;
open IN,$contig_file or die "Cannot open file $contig_file: $!\n";
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	$c[2]--;
	push @contig, [@c];
}
close IN;
# step 3
# read postions of intervals (CNV OR CNVR) to be tested.
my @test_len;
my %test;
open IN,$tested_file or die "Cannot open file $tested_file: $!\n";
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	if($tested_file =~ /\.bed/)
	{
		push @test_len, $c[2]-$c[1];
		
		$c[2]--;
	}
	else
	{
		push @test_len, $c[2]-$c[1]+1;
	}
	push @{$test{$c[0]}},[@c[1,2]];
}
close IN;
# step 4
# build empirical distribution of number of genomic features overlapping TESTED INTERVALS.
my @overlap_gf_dist;
for (1..$simu_num)
{
	open OUT,">$simu_out_prefix.$_";
	my %simu_test = ();
	my @simu_contig = @contig;
	my @simu_test_len = @test_len;
	
	%simu_test = simu(\@simu_test_len,\@simu_contig);

	for (keys %simu_test)
	{
		print OUT ">$_\n";
		for (@{$simu_test{$_}})
		{
			print OUT join("\t",@$_),"\n";
		}
	}
	close OUT;

	my $overlap_gf_num = cmp_pos(\%simu_test,\%gf);
	
	push @overlap_gf_dist,$overlap_gf_num;
}

# step 5
# write statistics on empirical distribution into file
open OUT,">$simu_out_prefix.stats";
my %overlap_gf_dist;
for (@overlap_gf_dist)
{
	$overlap_gf_dist{$_} ++;
}
print OUT "NumOverlappingGenomicFeatures\tFrequency\tPercentage\n";
my $dist_av;
for (sort {$a <=> $b} keys %overlap_gf_dist)
{
	print OUT $_,"\t",$overlap_gf_dist{$_},"\t";
	$overlap_gf_dist{$_} = $overlap_gf_dist{$_}/$simu_num;
	$dist_av += $_ * $overlap_gf_dist{$_};
	print OUT $overlap_gf_dist{$_}*100,"\n";
}

print "Average of the empirical distribution: $dist_av\n";
my $actual_overlap_gf_num = cmp_pos(\%test,\%gf);
print "Actual number of overlapping gf: $actual_overlap_gf_num\n";
my ($l_pvalue,$r_pvalue) = cal_pvalue(\%overlap_gf_dist,$actual_overlap_gf_num);
print "P-value: $l_pvalue $r_pvalue\n";


sub cal_pvalue {
	my ($ref_dist,$num) = @_;
	my $rp = 0;
	my $lp = 0;
	for (sort {$b <=> $a} keys %$ref_dist)
	{
		if($_ > $num)
		{
			$rp += $ref_dist->{$_};
		}
		elsif($_ < $num)
		{
			$lp += $ref_dist->{$_};
		}
	}
	
	return ($lp,$rp);
}

sub print_usage()
{
	print "Run the program:\n";
	print "  perl enrichment.pl  GenomicFeatures  contig.bed  WhatToBeTested.bed  Num_Replicates  Output_Prefix \n";
	print "  All files are tab-delimited file.\n";
	print "  GenomicFeatures: the first three columns are chrom, start and end\n";
	print "  contig.bed: a bed file containing contig positions\n";
	print "  WhatToTest.bed: a bed file containing positions of genomic sequences to test\n";
	print "  Num_Replicates: number of replicates to simulate\n";
	print "  Output_Prefix: name prefix of the output file\n";
	exit(1);
}

sub simu {
	my ($ref_len,$ref_contig)=@_;
	
	my %simu;
	
	while(@$ref_len)
	{
		my $rn1 = int rand @$ref_len;
		my $rand_len = $ref_len -> [$rn1];
	
		my @contig_seg = ();
		my $total_seg = 0;
		my $contig_num = 0;
		my $eff_contig_num = 0;
		my %eff2contig = ();
		for (@$ref_contig)
		{
			my @c = @$_;
			my $seg_num = $c[2]-$c[1]-$rand_len+2;
			if($seg_num > 0)
			{
				push @contig_seg, [$total_seg,$total_seg+$seg_num-1];
				$total_seg += $seg_num;
				$eff2contig{$eff_contig_num} = $contig_num;
				
				$eff_contig_num ++;
			}
			
			$contig_num ++;
		}
		my $rn2 = int rand $total_seg;
		$eff_contig_num = find_eff_contig_num($rn2,@contig_seg);
		
		$contig_num = $eff2contig{$eff_contig_num};
#		print $eff_contig_num,"\t",$contig_num,"\n\n";

		my ($contig_chrom,$contig_start,$contig_end) = @{$ref_contig->[$contig_num]};
		my ($contig_seg_start,$contig_seg_end) = @{$contig_seg[$eff_contig_num]};
		
		if($rn2 == $contig_seg_end)
		{
			push @{$simu{$contig_chrom}},[$contig_end-$rand_len+1,$contig_end];
			
			splice @$ref_contig,$contig_num,1,[$contig_chrom,$contig_start,$contig_end-$rand_len];
		}
		elsif($rn2 == $contig_seg_start)
		{
			push @{$simu{$contig_chrom}},[$contig_start,$contig_start+$rand_len-1];
			
			splice @$ref_contig,$contig_num,1,[$contig_chrom,$contig_start+$rand_len,$contig_end];
		}
		else
		{
			my $start_tmp = $contig_start+$rn2-$contig_seg_start;
			push @{$simu{$contig_chrom}},[$start_tmp,$start_tmp+$rand_len-1];
			
			splice @$ref_contig,$contig_num,1,[$contig_chrom,$contig_start,$start_tmp-1],[$contig_chrom,$start_tmp+$rand_len,$contig_end];
		}
		
		splice @$ref_len,$rn1,1;
	}
	
	return %simu;
}


sub cmp_pos {
	my ($ref_simu_test,$ref_gf)=@_;
	my $count = 0;
	for my $chrom (keys %$ref_gf)
	{
		next unless (defined $ref_simu_test->{$chrom});
		
		my @simu_chr_loc = @{$ref_simu_test->{$chrom}};
		@simu_chr_loc = sort {$a->[0] <=> $b->[0]} @simu_chr_loc;
		my @gf_chr_loc = @{$ref_gf->{$chrom}};
		@gf_chr_loc = sort {$a->[0] <=> $b->[0]} @gf_chr_loc;
		for (@gf_chr_loc)
		{
			my ($gf_s,$gf_e) = @$_;
			
			for (my $i=0;$i<@simu_chr_loc;$i++)
			{
				my ($simu_s,$simu_e) = @{$simu_chr_loc[$i]};
				
				if($simu_e < $gf_s)
				{
					splice @simu_chr_loc,$i,1;
					$i --;
				}
				elsif($simu_s > $gf_e)
				{
					last;
				}
				else
				{
					$count ++;
					last;
				}
			}
		}
	}
	return $count;
}

sub find_eff_contig_num {
	my ($pos,@test)=@_;
	my $j = $#test;
	my $i = 0;
	while(1)
	{
		my $x = int(($i+$j)/2);
		my $y = $x + 1;
		
		if($test[$x][0] <= $pos and $pos <= $test[$x][1])
		{
			return $x;
		}
		elsif($test[$y][0] <= $pos and $pos <= $test[$y][1])
		{
			return $y;
		}
		elsif($test[$x][0]>$pos)
		{
			$j = $x;
		}
		elsif($test[$y][1]<$pos)
		{
			$i = $y;
		}
	}
	
}
