# MPRAudit
Repository for MPRA analysis method "MPRAudit"

## Overview
MPRAudit calculates the fraction of variance explained by sequence variation, which we call b<sup>2</sup>.  In addition to code to simulate sample data, at its core MPRAudit consists of 5 functions to produce and return b<sup>2</sup>:

Estimate b<sup>2</sup>:

1. Calculate the technical variance

2. Calculate the total variance

Estimate var(b<sup>2</sup>):

3. Calculate the variance of the technical variance

4. Calculate the variance of the total variance

Then (5) combine the above results and return b<sup>2</sup> and var(b<sup>2</sup>).

There are minor variations depending on whether MPRAudit is being applied to the variation of individual sequences, pairs of sequences, or groups of sequences.

## Preparing Data
MPRAudit reads in data using pandas "read_csv" and assumes there is no header.

For individual sequences, the file should have 3 columns:

RNA_counts,DNA_counts,sequence_indicator

For groups of individual sequences, the file should have 4 columns:

RNA_counts,DNA_counts,sequence_indicator,group_indicator

For pairs of sequences, the file should have 6 columns:

RNA_counts1,DNA_counts1,RNA_counts2,DNA_counts2,sequence_indicator1,sequence_indicator2

Counts can be normalized (not necessarily integer).

"sequence_indicator" tells MPRAudit which reads are clones and which reads are distinct sequences.  For instance, if there are two sequences with three clones each, the sequence indicator column might be 1,1,1,2,2,2 or 15,15,15,14.9,14.9,14.9 (the actual values don't matter, MPRAudit looks for whether or not they're distinct or equal).

For pairs of sequences, "sequence_indicator" tells MPRAudit which sequences are paired, and they may have different numbers of clones.  For instance, there might be two pairs of sequences with different numbers of clones, and the data file might look like:\
6,7,6,7,1,1\
6,7,6,7,1,1\
6,6,6,8,1,1\
6,6,6,7,1,2\
6,7,6,7,2,2\
6,6,6,8,2,2\
7,7,,,2,


[NEED TO BE ABLE TO INCLUDE COLUMNS OF UNEQUAL LENGTHS]

By default, MPRAudit assumes the columns are tab-separated, but CSVs or other delimiters can be used by passing the -sep flag.

## Code
etc

## Citation
David A. Siegel, Olivier Le Tonqueze, Anne Biton, David J. Erle, and Noah Zaitlen, "MPRAudit Quantifies the Fraction of Variance Described by Unknown Features in Massively Parallel Reporter Assays" (2020).
