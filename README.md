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

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator

For groups of individual sequences, the file should have 4 columns: 

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator, (4) group_indicator

For pairs of sequences, the file should have 6 columns: 

(1) RNA_counts1, (2) DNA_counts1, (3) RNA_counts2, (4) DNA_counts2, (5) sequence_indicator1, (6) sequence_indicator2

Counts can be normalized (not necessarily integer).

"sequence_indicator" tells MPRAudit which reads are clones and which reads are distinct sequences.  For instance, if there are two sequences with three clones each, the sequence indicator column might be 1,1,1,2,2,2 or 15,15,15,14.9,14.9,14.9 (the actual values don't matter, MPRAudit looks for whether or not they're equal or distinct).

For pairs of sequences, "sequence_indicator" tells MPRAudit which sequences are paired, and they may have different numbers of clones.  For instance, there might be two pairs of sequences with different numbers of clones, and the data file might look like:\
6,7,6,7,1,1\
6,7,6,7,1,1\
6,6,6,9,1,1\
6,6,6,7,1,2\
6,7,6,7,2,2\
6,6,6,8,2,2\
7,7,,,2,\
6,6,,,2,

In this case the final rows have missing data because the number of clones differ.  The extra commas must be present for the data to load properly.  But there should be NO missing data in the middle of the file, e.g. this will lead to errors:\
6,7,6,7,1,1\
,7,6,7,1,1\
6,,6,9,1,1\
6,6,6,7,,2\
6,7,6,7,2,2\
6,6,6,8,2,2\
7,7,,,2,\
6,6,,,2,

This will cause errors because NaNs are removed after each column is loaded, so the program will mix up the clone counts like so:\
6,7,6,7,1,1\
6,7,6,7,1,1\
6,6,6,9,1,1\
6,7,6,7,2,2\
6,6,6,7,2,2\
7,7,6,8,2,2\
6,6,,,2,\


By default, MPRAudit assumes the files are CSVs (comma-separated), but other delimiters can be used by passing the -sep flag.  We don't recommend using white space (tabs/spaces) because columns may be unequal length and it might be beneficial to be able to observe missing data.

## Code
MPRAudit uses basic functionality from common python packages: pandas, numpy, and scipy.

MPRAudit functions are given in MPRAudit/MPRAudit_Functions.py
MPRAudit can be run as:

```
python MPRAudit.py -infile [input file] -outfile [output file] -ratio [ratio function indicator] -paired [True/False] -sep [delimiter] -timepoints [1 or 2] -numtrials [integer] -jackpow [float between 0 and 1]
```

The only required input is the input filename.
-If no output filename is given MPRAudit prints to standard output using the python "print" function.
-ratio uses RNA/(RNA+DNA) by default, must be an integer, currently 1-5 (see below for further details)
-paired is False by default
-sep is comma "," by default
-timepoints is 1 by default
-numtrials is 100 by default
-jackpow is 0.6 by default

## Citation
David A. Siegel, Olivier Le Tonqueze, Anne Biton, David J. Erle, and Noah Zaitlen, "MPRAudit Quantifies the Fraction of Variance Described by Unknown Features in Massively Parallel Reporter Assays" (2020).
