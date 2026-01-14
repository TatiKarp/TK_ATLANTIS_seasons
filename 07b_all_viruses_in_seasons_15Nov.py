#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 5 13:40 2024

This script will check the presence of different viruses in all atlantis samples 

@author: tatiana
"""
# %%
import csv
import os

#%%
# list all the files 
microb_dir = "/Users/tatiana/Work/RP2/ATLANTIS/Microbes/Output/nasal_brushes_score_check/"
files = os.listdir(microb_dir)

master_table = "/Users/tatiana/Work/RP2/ATLANTIS/Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"

samples_list = []
sample_season = {}
with open (master_table, 'r') as csvfile:                  
    reader = csv.DictReader(csvfile, delimiter =',') 
    for row in reader:
        samples_list.append(row['original_id'])
        sample_season[row['original_id']] = row['new_seasons']
        
# %%
# check viruses in a sample 
#file = "/Users/tatiana/Work/RP2/ATLANTIS/Microbes/Output/nasal_brushes_score_check/103893-001-373_z_score_check.csv"

all_taxa =  {} # all unique taxa IDs (both level 1 and 2) with z-scores >3 in a least 1 sample
all_samples = {}

# there is a problem with one sample 399 (I don't know what happened)
for sample in samples_list: 
    needed_file = "not found"
    for file in files:
        if sample in file:
            needed_file = microb_dir + file
    if needed_file != "not found":  
        with open (needed_file, 'r') as csvfile:         
            sample_name = sample
            print(sample_name)         
            reader = csv.DictReader(csvfile, delimiter =',') 
            all_samples[sample_name] = {} #per each sample create an empty dict to save the taxa later
            for row in reader:
                category = row['category']
                #print(category)
                if category == "viruses": 
                    taxa = row['taxa']
                    tax_level = row['tax_level']
                    genus_tax_id = row['genus_tax_id']
                    name = row['name']
                    is_phage = row['is_phage']
                    z_scores = []
                    for z_score in [row['new_nt_z_score'], row['new_nr_z_score']]:
                        try:
                            z_scores.append(float(z_score))
                        except ValueError:
                            pass
                    if z_scores: #if the list is not empty
                        max_z_score = max(z_scores) # get mx score from 
                    #save only if z_score >3 and it is not a phage
                        if max_z_score > 3 and is_phage == "false":
                            #print(z_scores, max_z_score)
                            all_samples[sample_name][taxa] = [tax_level, genus_tax_id, name, max_z_score]
                            all_taxa[taxa] = [tax_level, genus_tax_id, name]
    else:
        pass
            
            
# %% loop over taxa
# record N of samples per taxa that have the taxa present in 2 seasons

all_taxa_seasons = {}
for taxa in all_taxa:
    winter_spring = 0
    summer_autumn = 0
    all_taxa_seasons[taxa] = []
    for sample in all_samples:
        if taxa in all_samples[sample].keys():
            if sample_season[sample] == "winter_spring":
                winter_spring += 1
            elif sample_season[sample] == "summer_autumn":
                summer_autumn += 1
    all_taxa_seasons[taxa] = [winter_spring, summer_autumn]
    print(all_taxa_seasons[taxa])
            
# save only taxa that is present in at least 4 samples
at_least_4_samples_taxa_seasons = {}
n = 0
for taxa in all_taxa_seasons:
    if sum(all_taxa_seasons[taxa]) >=4:
        at_least_4_samples_taxa_seasons[taxa] = all_taxa_seasons[taxa]
        n += 1


#%% perform chi-suare test 
# winter_spring vs summer_autumn

from scipy.stats import chi2_contingency
import numpy as np
from collections import Counter

# create an array with n of samples in 2 seasons
group_counts = Counter(sample_season.values())
season_freq = np.array([group_counts.get("winter_spring", 0),group_counts.get("summer_autumn", 0)])


taxa_chi_test = {}
p_val_list = []

# loop over all the taxa
# data - contingency table seasons vs viral presence 
for taxa in all_taxa_seasons:
    data = [season_freq - np.array(all_taxa_seasons[taxa]),
            np.array(all_taxa_seasons[taxa])] 
    stat, p, dof, expected = chi2_contingency(data, correction = False)
    taxa_chi_test[taxa] = p
    p_val_list.append(p) # save all p-values for FDR correction
    if p < 0.05:
        print(data, taxa, p, all_taxa[taxa]) 

# %% perform FDR correction 
from statsmodels.stats.multitest import multipletests

p_values = np.array(p_val_list)
rejected, adjusted_p_values, _, _ = multipletests(p_values, method='fdr_bh')
min(adjusted_p_values)

# min fdr adjusted p-value == 0.466