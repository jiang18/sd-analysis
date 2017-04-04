# get_sets_of_similar_hits.pl
# This program is used to get sets of similar hits from the final output of 
# the WGAC pipeline.

use strict;
use warnings;
use 5.010;

@ARGV == 2 or &print_usage();

open IN,$ARGV[0] or die "Cannot open file $ARGV[0]: $!";    # the FINAL output of the WGAC pipeline

my @sd=();
# structure of @sd
# 
# $sd[] = {}
# Every element of @sd is a set of all similar hits, stored in a hash with keys being chromosomes.
#
# $sd[]->{chrom} = []
# Val in the hash is an array containing similar hits in the chromosome.
#
# $sd[]->{chrom}->[] = [start, end]
# The element in the array is further an array containing start and end of a hit.
#

say "Reading $ARGV[0]...";
# Collect all similar hits into one set
my $line_num;
MLOOP: while(<IN>)
{
	$line_num++;
	$line_num % 1000 == 0 and print "Finished processing $line_num lines.\n";
	chomp;
	my @c = split /\t/;
	$c[0] =~ s/.fa$//;
	$c[0] =~ s/^chrom//;
	
	$c[3] =~ s/.fa$//;
	$c[3] =~ s/^chrom//;
	
	# swap start and end of reverse alignments 
	$c[4] > $c[5] and @c[4,5]=@c[5,4];

	my $flag = 0;
	for my $sd_ref (@sd)
	{
		if(defined $sd_ref->{$c[0]})
		{
			for my $sd_loc (@{$sd_ref->{$c[0]}})
			{
				# determine if the current hit belongs to a set of similar hits
				# using approximate boundary.
				if(abs($c[1] - $sd_loc->[0])<1000 and abs($c[2] - $sd_loc->[1])<1000)
				{
					push @{$sd_ref->{$c[3]}},[@c[4,5]];
					$flag = 1;
					next MLOOP;
				}
			}
		}

		if(defined $sd_ref->{$c[3]})
		{
			for my $sd_loc (@{$sd_ref->{$c[3]}})
			{
				# determine if the current hit belongs to a set of similar hits
				# using approximate boundary.
				if(abs($c[4] - $sd_loc->[0])<1000 and abs($c[5] - $sd_loc->[1])<1000)
				{
					push @{$sd_ref->{$c[0]}},[@c[1,2]];
					$flag = 1;
					next MLOOP;
				}
			}
		}
	}
	
	unless($flag)
	{
		push @sd, {$c[0]=>[[@c[1,2]]]};
		push @{$sd[-1]->{$c[3]}},[@c[4,5]];
	}
}

say "Writing into output file";
open OUT,'>'.$ARGV[1];                              # Output
say OUT "NumChrom\tNumCopies\tHits(Chrom: Pos)";
# Count the number of distinct chromosomes and the number of copies for each of the sets
for my $sd_ref (@sd)
{
	my $chr_cnt = keys %$sd_ref; # number of distinct chromosomes
	my $copy_cnt; # number of copies
	for (keys %$sd_ref)
	{
		$copy_cnt += @{$sd_ref->{$_}};
	}

	print OUT $chr_cnt,"\t",$copy_cnt,"\t";
	for my $chrom (keys %$sd_ref)
	{
		print OUT $chrom,":";
		for(sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @{$sd_ref->{$chrom}})
		{
			print OUT " ";
			print OUT join "-",@$_;
		}
		print OUT "\t";
	}
	print OUT "\n";
}

sub print_usage()
{
	print "Run the program: \n";
	print "  perl get_sets_of_similar_hits.pl argv1 argv2 \n";
	print "  The first argument should be the final output file of the WGAC pipeline. \n";
	print "  The second argument should be the output file you want to specify. \n";
	exit;
}
