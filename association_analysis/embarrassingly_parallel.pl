# Please contact Jicai Jiang (jicai.jiang@gmail.com) if you have any questions or concerns.

use warnings;
use strict;
use POSIX ":sys_wait_h";

@ARGV==7 or die "Please run the program:  \n  perl embarrassingly_parallel.pl  parallel_start_number  parallel_end_number other_arguments.\n";

my ($parallel_start, $parallel_end, $gf_file, $contig_file, $tested_file, $simu_num, $simu_out_prefix)=@ARGV;
($parallel_start =~ /^\d+$/ and $parallel_end =~ /^\d+$/) or die "The first 2 argvs must be integer.\n";

my $zombies = 0;
my $num_proc = 0;
my $collect;
$SIG{CHLD} = sub { $zombies++; $num_proc-- };

foreach my $i ($parallel_start .. $parallel_end)
{
	my $pid = fork();
	if(!defined($pid))
	{
		print "Error in fork $i: $!";
		exit 1;
	}
	
	if($pid == 0)
	{
		print "Child $i (PID = $$): start\n";
		system("perl enrichment.pl $gf_file $contig_file $tested_file $simu_num $i.$simu_out_prefix");
# Finished
		print "Child $i (PID = $$): end\n";
		exit(0);
	}
	$num_proc++;
	if($zombies > 0)
	{
		while (($collect = waitpid(-1, WNOHANG)) > 0)
		{
			$zombies --;
		}
	}
	while($num_proc >= 40)
	{
		sleep(1);
	}
}
exit 0;
