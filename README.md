# Cryptococcus_more_FPs_2025

This code is to generate more codon-adjusted fluorescent proteins for Cryptococcus neoformans.

Edward Wallace, Edward.Wallace@ed.ac.uk, 2025-03-21.

Funded by the Wellcome Trust Bioimaging award.

(Note: README written post-hoc in 2026.)

Related repositories generating FPs for other fungal species:

- https://github.com/ewallace/Yeast_FPs_2025
- https://github.com/ewallace/Emergomyces_constructs_2025

# Contents

## data

Input data:

- codon frequency table for C. deneoformans JEC21 5% highest-translated genes, `codFreqTableJEC21HiTrans.txt`. From [Wallace, Maufrais et al. 2020, doi:10.1093/nar/gkaa060](http://dx.doi.org/10.1093/nar/gkaa060)
- protein sequences for fluorescent proteins from fpbase.org. Downloaded by EW February 2025.

Red FPs:

- [mCherry-XL](https://www.fpbase.org/protein/mcherry-xl/)
- [mKO2, aka mKusabira-Orange2](https://www.fpbase.org/protein/mko2/)
- [mNeptune2.5](https://www.fpbase.org/protein/mneptune25/), but the filenames are called mNeptune2-5 because filenames don't like dots
- [mScarlet-I3](https://www.fpbase.org/protein/mscarlet-i3/)
- [mScarlet3-H, aka mYongHong](https://www.fpbase.org/protein/mscarlet3-h/)


Green & yellow FPs:

- [mGreenLantern](https://www.fpbase.org/protein/mgreenlantern/)
- [fuGFP, aka Free Use GFP](https://www.fpbase.org/protein/free-use-gfp/)
- [hfYFP, aka Hyperfolder YFP](https://www.fpbase.org/protein/hyperfolder-yfp/)
- [mStayGold, aka QC2-6 FIQ](https://www.fpbase.org/protein/mstaygold/)

HaloTag fluorogenic fluroescent protein:

- HaloTag, sequence checked against [YIplac211-Sec31-yoHalo
(Addgene Plasmid #115423)](https://www.addgene.org/115423/) and snapgene annotation of HaloTag

Linker sequences:

- GGGGSGS, a flexible amino acid linker


## designs

Coding sequences for fluorescent proteins, adjusted to 

- codon frequencies of Histoplasma capsulatum ribosomal proteins, 
- moderate GC content, and 
- to avoid select restriction enzyme sites. 

See source code for details.

## src

Source code in python format. Uses the DNAchisel package and Biopython tools.

- codon_adjust_FPs.py - using DNAchisel's DnaOptimizationProblem to design sequences encoding fluorescent proteins
- codon_count_functions.py - functions to calculate codon counts and frequencies


# How to run the code

To runs the code to design the proteins, from the root directory of of `Cryptococcus_more_FPs_2025`, in a bash terminal run:

```bash
python src/codon_adjust_FPs.py
```

## Computational environment

You'll need the relevant python packages installed. Biopython and DNAchisel, etc.

It should be possible to install the packages using conda:

```bash
conda create --name codons2025 --file codons2025_env.txt
conda activate codons2025
```

Beware, I have not tested this conda installation.
