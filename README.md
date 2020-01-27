# MPRAudit
Repository for MPRA analysis method "MPRAudit"

## Code To Generate Figures in the Paper
I've posted code to generate the figures in the paper, which you can use as examples for how to run MPRAudit.  In addition to code to simulate sample data, MPRAudit consists of 5 simple functions to produce and return b<sup>2</sup>:

## Functions:
### Function "return_counts"
Simulates RNA and DNA counts for a given number of clones.  Takes as input several experimental parameters (number of clones, amount of noise, etc).  Returns as output a list of RNA counts for each clone, a list of DNA counts for each clone, and variables proportional to the total number of counts of RNA and DNA (which may be useful for normalization).

### Function "return_counts_T4T0"
Simulates RNA and DNA counts for a given number of clones at two time points.


### Function "deleteDjackknife_variance_T0"
Takes as input two array-like objects that represents counts of RNA and DNA (respectively), returns 

