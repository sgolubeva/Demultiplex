Demultiplexing results

General statistics results 

|      |Count|   Percent|
|------|------|----------|
|Total|   363246735|        100.00|
|Matched| 331755033|        91.33|
Index Hopped|    707740|   0.19|
Unknown| 30783962|         8.47|

Overall amount of index swapping in my results was 0.19% which is a relatively small number.
Illumina claims that levels of index swapping around 1.5% is acceptanle for their patterned 
flowcells. There were two indexes that showed higher number of hopping events than others. I can attribute
that to higher concentration of free floating adapters that were not cleaned up during the library prep.

The detailed metrics for each index combination index hopping events could be found in index_hopped_stats.tsv file

Additionally, there are counts and percents for matched indexes that could be found in matched_percents.tsv

Overall, if the library pooling was uniform it would have resulted in close to equal counts for each matched index.
The counts I obtained are not uniform with some indexes getting more counts than others which I could attribute to 
poor library pooling.

For quality score filtering, I only filter the indexes that have 'N's. I am not implementing any q-score cutoff at this time
for reads because filtering Ns should be sufficient for our purposes. In the future, I want to implement a q-score cutoff defined by the user as an optional argument in the argparse. For our data, the smallest score is 12, which is a little higher than 90%
correct. It will take 3 wrong calls to change one index into another index from our list. Also the changes should happen to both indexes to be converted into another index from the list. The likelyhood of that is very small, that's why I think that filtering only Ns is sufficient. 