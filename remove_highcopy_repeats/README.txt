# Remove high-copy repeats based on the following criteria:
# 1) >= 50 copies in the genome
# OR
# 2) a distribution to multiple chromosomes (>= 3)

# Required file: the final output of the WGAC pipeline

Here is an example on how to use the programs.
# Three steps
# 1
perl get_sets_of_similar_hits.pl ./test_data/test.wgac ./test_data/test.raw_hit_sets
# 2
perl clean_sets_of_similar_hits.pl ./test_data/test.raw_hit_sets ./test_data/test.clean_hit_sets
# 3
perl rm_highcopy_repeats.pl ./test_data/test.clean_hit_sets ./test_data/test.wgac ./test_data/test.wgac.clean

# Final output: ./test_data/test.wgac.clean

# Note that the programs may generate conservatively cleaned results, since 
# the number of copies may be overestimated. Generally, this does not matter, 
# because the overestimation of copy numbers can be often seen for high-copy 
# repeats but rarely seen for low-copy ones.

# Please contact Jicai Jiang (jicai.jiang@gmail.com) if you have any questions or concerns.
