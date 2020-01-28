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

This must be tweaked depending on whether MPRAudit is being applied to single sequence fragments, or pairs or groups of sequences.

## Preparing Data
MPRAudit expects data of the form [describe, then show vectors or something]

## Code
etc

## Citation
