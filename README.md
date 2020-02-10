# MPRAudit
Repository for MPRA analysis method "MPRAudit"

## Overview
MPRAudit takes a CSV as input; it returns the fraction of variance explained by sequence variation, which we call b<sup>2</sup>, as well as its estimated uncertainty, the estimated standard deviation of b<sup>2</sup>.  At its core MPRAudit consists of 5 functions:

Estimate b<sup>2</sup>:

1. Calculate the technical variance

2. Calculate the total variance

Estimate var(b<sup>2</sup>):

3. Calculate the variance of the technical variance

4. Calculate the variance of the total variance

Then (5) combine the above results to return b<sup>2</sup> and SD(b<sup>2</sup>).  We return the standard deviation rather than the variance because this is what most people care about in practice.

There are minor variations depending on whether MPRAudit is being applied to the variation of individual sequences, pairs of sequences, or groups of sequences.

We also provide code to generate sample data, example scripts, example data, and the python notebook that generated the figures in our manuscript.

## Types of Data

A. MPRAs like Fast-UTR (Zhou et al, Nature Biotech 2014) return counts of RNA and DNA.

B. Some MPRAs might be interested in the amount of variation remaining after groups of sequences are taken into account.

C. Some MPRAs like our recent bioRxiv submission "Massively Prallel Analysis of Human 3' UTRs Reveals that AU-Rich Element Length and Registration Predict mRNA Destabilization" might measure counts of RNA and DNA at separate time points, and might be interested in a ratio of ratios, such as RNA/(RNA+DNA) evaluated at T4/(T4+T0).

D. Some MPRAs such as our recent bioRxiv submission might be interested in pairs of reference and mutant sequence data (how much of the variation in the data is due to the mutations as opposed to technical variance?).

E. MPRAs can be at multiple time points and have pairs of reference and mutant sequences.

F. Some assays like CRISPR screens return only one type of data, such as RNA counts.

## Preparing Data
MPRAudit reads in data using pandas "read_csv" and assumes there is no header.  Note that counts do not need to be integer, they can be normalized or transformed.

A. For standard MPRA data with RNA counts and DNA counts, the input file should have 3 columns and no special flags are required: 

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator

B. For calculating b<sup>2</sup> within groups of sequences, the file should have 4 columns and no special flags are required: 

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator, (4) group_indicator

C. For MPRA data at two time points, the file should have 6 columns and "-timepoints 2" should be used:

(1) RNA_counts_T0, (2) DNA_counts_T0, (3) RNA_counts_T4, (4) DNA_counts_T4, (5) sequence_indicator_T0, (6) sequence_indicator_T4

D. For pairs of MPRA sequence data, the file should have 6 columns and "-paired True" should be used:

(1) RNA_counts1, (2) DNA_counts1, (3) RNA_counts2, (4) DNA_counts2, (5) sequence_indicator1, (6) sequence_indicator2

For pairs of sequences as in (C) or (D), "sequence_indicator" tells MPRAudit which sequences are paired, and they may have different numbers of clones.  For instance, there might be two pairs of sequences with different numbers of clones, and the data file might look like:\
6,7,6,7,1,1\
2,3,6,8,1,1\
6,7,6,5,1,1\
3,4,6,10,1,2\
6,7,6,15,2,2\
16,10,6,10,2,2\
7,7,,,2,\
6,6,,,2,

In this case, there are two sets of RNA and DNA data, with two sequences labelled 1 and 2 given in the 5th and 6th columns.  For the first set of RNA and DNA data, sequence 1 has 4 clones and sequence 2 has 4 clones (there are 4 1's and 4 2's in the 5th column).  For the second set of RNA and DNA data, sequence 1 has 3 clones and sequence 2 has 3 clones (there are 3 1's and 3 2's in the 6th column).

We saved this small dataset as "ExampleData1.csv" in MPRAudit/Examples/.  To run MPRAudit on this as paired RNA/(RNA+DNA) data, the command is
```
python MPRAudit.py -infile /path/to/ExampleData1.csv -paired True
```
and the output is roughly (results vary due to randomization; use -numtrials for more precise values):
```
b2_mean: 0.83
b2_std: 0.07
```

The paired b<sup>2</sup> is large because the clones of sequence 1 are comparable in size for both pairs, but the fourth column is much larger than the third column for sequence 2 while the second column is smaller on average than the first column for sequence 2.


In this case the final rows have missing data because the number of clones differ.  The extra commas must be present for the data to load properly.  But there should be NO missing data in the middle of the file, e.g. this will lead to errors and an exception might not be thrown:\
6,7,6,7,1,1\
2,3,6,8,1,1\
6,7,6,5,1,1\
3,4,,,1,\
6,7,6,10,2,2\
16,10,6,15,2,2\
7,7,6,10,2,2\
6,6,,,2,


E. For pairs of data at two time points, the file should have 12 columns: (1) RNA_counts1_T0, (2) DNA_counts1_T0, (3) RNA_counts1_T4, (4) DNA_counts1_T4, (5) RNA_counts2_T0, (6) DNA_counts2_T0, (7) RNA_counts2_T4, (8) DNA_counts2_T4, (9) sequence_indicators1_T0, (10) sequence_indicators1_T4, (11) sequence_indicators2_T0, (12) sequence_indicators2_T4


F. Finally, for CRISPR screens and other assays where only RNA counts are obtained and no DNA counts, the input file should have 2 columns and you must use "-CRISPR_log_flag True/False" to tell MPRAudit whether to compare the raw values or log2(counts+1):

(1) RNA_counts, (2) sequence_indicator

By default, MPRAudit assumes the files are CSVs (comma-separated), but other delimiters can be used by passing the -sep flag.  We don't recommend using white space (tabs/spaces) because columns may be of unequal length and it might be beneficial to be able to observe missing data.

## Code
MPRAudit uses basic functionality from common python packages: pandas, numpy, and scipy.

MPRAudit functions are given in MPRAudit/MPRAudit_Functions.py.  To make sure MPRAudit runs properly, make sure MPRAudit_Functions.py and MPRAudit.py are in the same folder.

MPRAudit can be run as:

```
python MPRAudit.py -infile [input file] -outfile [standard out] -ratio [1] -paired [False] -sep [","] -timepoints [1] -numtrials [100] -jackpow [0.6] -CRISPR_log_flag [None]
```

* -infile requires the input file (and path if appropriate)
* -outfile gives MPRAudit an output filename (it prints to standard output by default)
* -ratio uses RNA/(RNA+DNA) by default, must be an integer, currently 1-5 for paired RNA/DNA data, 11-13 for two time points (see below for further details)
* -paired is False by default, whether or not paired differences between sequences are being examined
* -sep is comma "," by default, the text delimiter in the input file (commas in CSVs, "\t" in tab-delimited files, etc)
* -timepoints is 1 by default, option for 2 exists
* -numtrials is 100 by default, more gives better output accuracy
* -jackpow is a value between 0 and 1 and tells the delete-D jackknife what proportion of sequences to hold out
* -CRISPR_log_flag is None by default, and should only be used if only one type of data is being analyzed (like RNA data in a CRISPR screen as opposed to RNA/DNA data in an MPRA)

### Ratio Functions
There are currently several choices of ratio functions.  For RNA and DNA data:
1. RNA/(RNA+DNA)
2. log2(RNA/DNA)
3. log2((RNA+1)/(DNA+1))
4. RNA/DNA
5. (RNA+1)/(DNA+1)

For RNA and DNA data at two time points (T4 and T0 in our publication):

11. T4/(T4+T0) where T0 = RNA/(RNA+DNA)|T0 and T4 = RNA/(RNA+DNA)|T4\
12. T4/T0 where T0 = log2(RNA/DNA)|T0 and T4 = log2(RNA/DNA)|T4\
13. T4/T0 where T0 = log2(RNA/DNA+1)|T0 and T4 = log2(RNA/DNA+1)|T4


Feel free to implement your own and/or contact us with suggestions!


### Additional Example Simulations

We provided a small CSV as discussed above to test the basic function of MPRAudit on paired data with a CSV input.

To test MPRAudit without creating a CSV file, some additional example simulations can be run as

```
python Example_Simulation#.py
```

On a 2014 Mac Mini,\
Simulation1 returns 7 seconds, b2 ~ 0.51, SD ~ 0.03\
Simulation2 returns 14 seconds, b2 ~ 0.88, SD ~ 0.01\
Simulation3 returns 13 seconds, b2 ~ 0.35, SD ~ 0.04\
Simulation4 returns 6 seconds, b2 ~ 0.53, SD ~ 0.03


## Citation
David A. Siegel, Olivier Le Tonqueze, Anne Biton, David J. Erle, and Noah Zaitlen, "MPRAudit Quantifies the Fraction of Variance Described by Unknown Features in Massively Parallel Reporter Assays" (2020).

## Contact
David [dot] Siegel [at] ucsf [d0t] edu, please put "MPRAudit" in the subject line.  This is a work in progress.

## Thanks
DAS would like to thank Christa Caggiano for help with this Github
