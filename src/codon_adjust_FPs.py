#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2025-03-21

@author: Edward Wallace
"""
# This script uses DNAchisel to produces Saccharomyces cerevisiae codon-adjusted DNA sequences from an amino acid sqeuence
# Instruction for DNAchisel can be found at https://edinburgh-genome-foundry.github.io/DnaChisel/
# DNAchisel paper: https://doi.org/10.1093/bioinformatics/btaa558
# This script builds on anidulans_codon_adjust.py script from Domenico Modaffari


# import codon counting functions, EW adapted from CodonAdaptationIndex package
from codon_count_functions import *

# import relevant functions from DNAchisel
from dnachisel import (DnaOptimizationProblem,
                       EnforceTranslation,
                       CodonOptimize,
                       reverse_translate,
                       EnforceGCContent,
                       UniquifyAllKmers,
                       AvoidRareCodons,
                       AvoidPattern, #use this if you want to exclude restriction enzymes cut sites from DNA sequence
                       AvoidHairpins
                       )

# import SeqFeature, for labeling
from Bio import SeqFeature as sf

# import date, so we can time-stamp the output
from datetime import date

# Ensure reproducible output by setting the seed for the random number generator
# This is because the DnaOptimizationProblem uses some random processes in its optimization
from numpy import random

# Import pandas for data frames
import pandas as pd

# Calculate the codon frequencies to use
CdCodonFreqTable = pd.read_table("data/codFreqTableJEC21HiTrans.txt")

CdCodonFreqs = CodonFrequencyByAA(
    counts = pd.Series(CdCodonFreqTable.Count.values, index=CdCodonFreqTable.Codon).to_dict(),
        outputtype = "freq"
    )
    

# Define constraints to avoid
# restriction enzyme sites, that rely on Bio.Restriction as called in dnachisel.AviodPattern
enzymes_to_avoid = [ 'BamHI', 'NdeI', 'HindIII', 'SwaI', 'AgeI', 'XhoI', 'SpeI', 'KpnI', 'NotI', 'EcoRI', 'BglII', 'SacI', 'SacII', 'BsaI', 'BsmBI', 'AarI']
sites_to_avoid = [ AvoidPattern("%s_site" % enzyme)
    for enzyme in  enzymes_to_avoid ]

# also ensure translation, medium GC content, avoid rarest codons, avoid hairpins
other_constraints = [
    EnforceTranslation(),
    EnforceGCContent(mini=0.3, maxi=0.7, window=40), #GC content setting
    AvoidRareCodons(min_frequency = 0.05, codon_usage_table = CdCodonFreqs), #do not use codons below frequency of 0.1
    AvoidHairpins(stem_size=15) #avoid hairpins
    ] 

all_constraints = sites_to_avoid + other_constraints

# define objectives
objectives_match_uniquify = [
    CodonOptimize(method = "match_codon_usage", codon_usage_table = CdCodonFreqs), #adjust codons
    UniquifyAllKmers(k = 7) #try to avoid repeats in sequence, makes gene synthesys easier
    ]

# function that does all the adjustment/optimisation for a protein, and writes out to a record
def codon_adjust(proteinseq, 
        seed = None, 
        filepath = None, 
        record_id = None, 
        CDS_id = None,
        add_date = True,
        with_constraints = False,
        with_objectives = False,
        max_iters = 1000):
    
    # set random seed before starting the optimization problem
    if (seed is not None) :
        random.seed(seed)
    
    # This is the heart of the script
    problem = DnaOptimizationProblem(
        sequence = reverse_translate(proteinseq),
        constraints = all_constraints,
        objectives = objectives_match_uniquify
        )
    
    problem.max_random_iters = max_iters
    problem.resolve_constraints()
    problem.optimize()
    
    # After this, it's about labeling/annotating the output
    # convert to record
    record = problem.to_record(
        record_id = record_id,
        with_constraints = with_constraints,
        with_objectives = with_objectives)
    
    if (CDS_id is not None):
        # add the CDS feature into the record
        CDSfeature = sf.SeqFeature(
           location = sf.FeatureLocation(0, len(record.seq), strand = +1), 
           type = "CDS",
           id = CDS_id,
           qualifiers = {'gene':CDS_id})
        record.features = record.features + [ CDSfeature ] 
    
    if (add_date):
        # add date to record
        # note this doesn't write the date to the output file, don't know why
        record.annotations["date"] =  str(date.today())
    
    if (filepath is not None) :
        SeqIO.write(record, handle = filepath, format = "genbank")
    return record




## Apply to the proteins of interest. Write out to a record
# example if we were doing for single protein:
# import pdb; pdb.set_trace();
# codon_adjust(
#     proteinseq = SeqIO.read("data/mCardinal.fasta", "fasta"),
#     seed = 1,
#     filepath = "designs/mCardinaltest_2025-03-21.gbk",
#     record_id = "mCardinaltest",
#     CDS_id = "mCardinal",
#     max_iters = 10000)

# Loop over the relevant proteins:
proteins_for_design = [
    # "mCardinal",
    "mCherry-XL",
    # "mCherry",
    "mNeptune2-5",
    "mScarlet-I3",
    "mScarlet3-H",
    "mGreenLantern",
    "mKO2",
    # "mNeonGreen",
    "fuGFP",
    "hfYFP",
    "mStayGold",
    "GGGGSGS",
    "HaloTag"
    ]

for protein in proteins_for_design :
    print("codon adjusting for protein " + protein)
    # make a random seed for each protein based on the name, so reproducible
    pseed = int.from_bytes(protein.encode()) % 1000000
    # do the codon adjusting and write to file:
    codon_adjust(
        proteinseq = SeqIO.read("data/" + protein + ".fasta", "fasta"),
        seed = pseed,
        filepath = "designs/" + protein + "-Cno_2025-03-21.gbk",
        record_id = protein + "-Cno",
        CDS_id = protein)

print("all done")
