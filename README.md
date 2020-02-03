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

Counts can be normalized (not necessarily integer).
"sequence_indicator"

## Code
etc

## Citation
David A. Siegel, Olivier Le Tonqueze, Anne Biton, David J. Erle, and Noah Zaitlen, "MPRAudit Quantifies the Fraction of Variance Described by Unknown Features in Massively Parallel Reporter Assays" (2020).
