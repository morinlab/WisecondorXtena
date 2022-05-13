# Improved_WisecondorX
A modified version of the copy number caller, WisecondorX (https://github.com/CenterForMedicalGeneticsGhent/WisecondorX), that now infers purity and produces improved copy number calls from lpWGS data.

## WisecondorXtena :wink:
This version of WisecondorX utilizes the original code with no specifications other than bin size, calculates purity and how offset the mean log2 ratio is from zero, and reiterates these values back into the modified code as parameters that adjusts the genomic segments accordingly.

### The original WisecondorX is consistent of 3 main steps
* File conversion (aligned reads (BAM) to NPZ file format)
* Reference creation (using reference NPZ files)
* CNV prediction (through a user-definable cut-off ("beta"), or the default (leverages Z-scores))


### The modified WisecondorX is run via snakefile, consisting of two sections
**The run_updated_wisecondorX.smk set up:** 

Section 1 (initial run to calculate offset mean and purity):
* Runs the 3 WisecondorX steps with the original parameters and specifies 500kb bins.
* Plots the copy number profile that would be originally output by WisecondorX prior to any modifications made (this is called "unshifted") along with the corresponding density plot from a Gaussian mixture model.
* Determines how offset the mean log2 ratio is from zero (corresponding to how biased the assigned genomic aberrations are from the copy neutral state), as well as estimates the purity of each sample.

Section 2 (second run of WisecondorX to adjust segments):
* Reruns WisecondorX the offset value and purity estimate and inputs those as the offset_neut_peak and beta paramters, respectively.
* Produces new copy number profile plots with the density curves, after the adjustment is performed.
* A new purity estimate is calculated. This can be compared to IchorCNAs tumour fraction value.
  * Note the offset line in the density curve and the offset mean value are both no longer relevant after shifting is complete and the second set of plots are generated.

## New parameters (See the WisecondorX github page for the pre-existing paramaters and their description)
**offset_neut_peak** - A value calculated by the mixture model (using 1 cluster) that can be positive or negative (positive will cause an upwards shift and negative will cause a downwards shift in genomic segments). 
* How it works : this value is the mean log2 ratio of the neutral peak from the mixture model. It is added to the results_r, results_z and results_w variables to shift resulting log2 ratios and plots. 
