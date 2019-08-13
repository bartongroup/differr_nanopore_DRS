## A tool for detecting modifications from Nanopore DRS errors using a low modification control

### How to use:
* Basecall & map direct RNA data (all datasets **MUST** be sequenced at approximately the same time with the exactly the same flowcell/kit/MinKNOW software, and basecalled and mapped in the same way with the same model, otherwise you **WILL** get loads of false positives. Nanopore update their pores/kits/models/software all the time, so be wary...).
* Run `differr` on mapped bam files.
* Output from `differr` is a bed file with the positions with significantly altered error, and a optional hdf5 file with all of the per reference base basecalls, which might be useful for further analyses.
* Columns of bed file are:
    * chrom, start, end, name
    * score: -log10 of the FDR, rounded to nearest whole number
    * strand
    * odds ratio: the change in the ratio of matches to mismatches in the wild type compared to the mutant with low modifications. An odds ratio > 0 indicates more modifications in the WT.
    * G statistic for the comparison of pooled WT and mutant samples.
    * -log10 P value for the comparison of pooled WT and mutant samples.
    * -log10 FDR for the comparison of pooled WT and mutant samples.
    * G statistic for the homogeneity test of mutant replicates.
    * -log10 P value for the homogeneity test of mutant replicates.
    * G statistic for the homogeneity test of wild type replicates.
    * -log10 P value for the homogeneity test of wild type replicates.