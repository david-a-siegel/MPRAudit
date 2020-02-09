# MPRAudit
Repository for MPRA analysis method "MPRAudit"

## Overview
MPRAudit calculates the fraction of variance explained by sequence variation, which we call b<sup>2</sup>.  At its core MPRAudit consists of 5 functions to produce and return b<sup>2</sup>:

Estimate b<sup>2</sup>:

1. Calculate the technical variance

2. Calculate the total variance

Estimate var(b<sup>2</sup>):

3. Calculate the variance of the technical variance

4. Calculate the variance of the total variance

Then (5) combine the above results to return b<sup>2</sup> and SD(b<sup>2</sup>).  We return the standard deviation rather than the variance because this is what most people care about in practice.

There are minor variations depending on whether MPRAudit is being applied to the variation of individual sequences, pairs of sequences, or groups of sequences.

We also provide code to generate sample data, example scripts, and the python notebook that generated the figures in our manuscript.


## Preparing Data
MPRAudit reads in data using pandas "read_csv" and assumes there is no header.

For individual sequences, the file should have 3 columns: 

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator

For groups of individual sequences, the file should have 4 columns: 

(1) RNA_counts, (2) DNA_counts, (3) sequence_indicator, (4) group_indicator

For pairs of sequences, the file should have 6 columns: 

(1) RNA_counts1, (2) DNA_counts1, (3) RNA_counts2, (4) DNA_counts2, (5) sequence_indicator1, (6) sequence_indicator2

Counts can be normalized (not necessarily integer).

For pairs of sequences, "sequence_indicator" tells MPRAudit which sequences are paired, and they may have different numbers of clones.  For instance, there might be two pairs of sequences with different numbers of clones, and the data file might look like:\
6,7,6,7,1,1\
2,3,6,8,1,1\
6,7,6,5,1,1\
3,4,6,10,1,2\
6,7,6,15,2,2\
16,10,6,10,2,2\
7,7,,,2,\
6,6,,,2,

In this case, there are two sets of RNA and DNA data, with two sequences labelled 1 and 2 given in the 5th and 6th columns.  For the first set of RNA and DNA data, sequence 1 has 4 clones and sequence 2 has 4 clones (there are 4 1's and 4 2's in the 5th column).  For the second set of RNA and DNA data, sequence 1 has 3 clones and sequence 2 has 3 clones (there are 3 1's and 3 2's in the 6th column).

We saved this small dataset as "ExampleData1.csv".  To run MPRAudit on this as paired RNA/(RNA+DNA) data, the command is
```
python MPRAudit.py -infile /path/to/ExampleData1.csv -paired True
```
and the output is roughly (results vary due to randomization):
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




By default, MPRAudit assumes the files are CSVs (comma-separated), but other delimiters can be used by passing the -sep flag.  We don't recommend using white space (tabs/spaces) because columns may be of unequal length and it might be beneficial to be able to observe missing data.

## Code
MPRAudit uses basic functionality from common python packages: pandas, numpy, and scipy.

MPRAudit functions are given in MPRAudit/MPRAudit_Functions.py

MPRAudit can be run as:

```
python MPRAudit.py -infile [input file] -outfile [output file] -ratio [ratio function indicator] -paired [True/False] -sep [delimiter] -timepoints [1 or 2] -numtrials [integer] -jackpow [float between 0 and 1]
```

* -infile requires the input file (and path if appropriate)
* -outfile gives MPRAudit an output filename (it prints to standard output by default)
* -ratio uses RNA/(RNA+DNA) by default, must be an integer, currently 1-5 (see below for further details)
* -paired is False by default
* -sep is comma "," by default
* -timepoints is 1 by default
* -numtrials is 100 by default
* -jackpow is 0.6 by default

### Ratio Functions
There are currently five ratio functions:
1. RNA/(RNA+DNA)
2. log2(RNA/DNA)
3. log2((RNA+1)/(DNA+1))
4. RNA/DNA
5. (RNA+1)/(DNA+1)

Feel free to implement your own by editing "ratio_function" in MPRAudit_Functions and/or contact us with suggestions!


### Example Simulations

Example simulations can be run as

```
python Example_Simulation#.py
```

On a 2014 Mac Mini,\
Simulation1 returns 7 seconds, b2 ~ 0.51, SD ~ 0.03\
Simulation2 returns 14 seconds, b2 ~ 0.88, SD ~ 0.01\
Simulation3 returns 13 seconds, b2 ~ 0.35, SD ~ 0.04

## Citation
David A. Siegel, Olivier Le Tonqueze, Anne Biton, David J. Erle, and Noah Zaitlen, "MPRAudit Quantifies the Fraction of Variance Described by Unknown Features in Massively Parallel Reporter Assays" (2020).

## Contact
David [dot] Siegel [at] ucsf [d0t] edu, please put "MPRAudit" in the subject line

## Thanks
DAS would like to thank Christa Caggiano for help with this Github
