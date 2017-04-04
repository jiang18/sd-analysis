# Test if specific genomic sequences (like CNVs) are enriched in genomic features 
# (like gene regions, SDS) through simulation.

# The command to analyze the data in the example folder:
perl enrichment.pl ./example/genomic_features.txt ./example/contig.bed ./example/tested.bed 10 ./example/test

# To make use of multiple processors, use the script for embarrassingly parallel:
perl embarrassingly_parallel.pl 1 20 ./example/genomic_features.txt ./example/contig.bed ./example/tested.bed 10 ./example/test

# Please contact Jicai Jiang (jicai.jiang@gmail.com) if you have any questions or concerns.
