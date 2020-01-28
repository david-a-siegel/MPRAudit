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

5. Take the above outputs and return b<sup>2</sup> and var(b<sup>2</sup>).

## Test



1. "deleteDjackknife_variance_T0" returns the technical variance
2. np.var returns the total variance
3. "ordinary_jackknife_variance" returns the variance of the technical variance
4. "ordinary_variance_jackknife_variance" returns the variance of the total variance
5. "uncertainty" takes the above outputs and returns b<sup>2</sup> and var(b<sup>2</sup>)

There are also functions that apply to groups and pairs of sequence fragments or perturbations

## Functions:
### Function "return_counts"
Simulates RNA and DNA counts for a given number of clones.  Takes as input several experimental parameters (number of clones, amount of noise, etc).  Returns as output a list of RNA counts for each clone, a list of DNA counts for each clone, and variables proportional to the total number of counts of RNA and DNA (which may be useful for normalization).

### Function "return_counts_T4T0"
Simulates RNA and DNA counts for a given number of clones at two time points.


### Function "deleteDjackknife_variance_T0"
Takes as input two array-like objects that represents counts of RNA and DNA (respectively), returns 

