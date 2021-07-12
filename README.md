# SimpleScripts
simple scripts I wrote for one-off uses.

<b>HammingSimilarityMatrix.py</b>
This performs an all-vs-all comparison for same-length aligned sequences in multifasta format. 
Returns a matrix with a similarity metric based on Hamming distance, along the lines of "percent identity"

<b>nucleotide_diversity_metrics.py</b>
Calculates two different metrics for diversity for an input multifasta, again assumed to be aligned sequences of same length. 
First metric is "pi" - or the average percentage of positions that are different.
Second metric is average positions different per 100nt of total length.
Gaps (aka indels) are treated as differences, just the same as substitutions.

<b>ATassignment.py</b>
Parses an input aligned multifasta and returns only nonredundant sequences with numbered "allelic type" assignments.
