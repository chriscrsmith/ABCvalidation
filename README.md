# ABCvalidation
Code used for "Leveraging whole genome sequencing data for demographic inference with approximate Bayesian computation" (Smith and Flaxman, Mol Ecol Resources, 2019).

msABC_haplotypeStats contains the modified source code for msABC, plus the haplotype statistics described in our publication. See Makefile.

msp_wrapper.py and slim_wrapper.py are for building msprime and SLiM commands implementing user-specified prior ranges. 

slim2ms.py converts SLiM output to ms format.

SLiMmodels/ contains the edos code for the slim models used in our publication.

The "comparePhased" R scripts are what I used to analyze msABC output.





