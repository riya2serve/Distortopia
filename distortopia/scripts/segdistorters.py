"""Identify variants showing distorted segregation patterns."""
import pandas as pd

def detect(args, csv = False):
    genotypes = parse_vcf(args.f1)
    exp_ratio = 0.5 #expected ratio
    distorted = []

    if csv:
        obs_exp = {"obs": [], "exp": []} #Adding an empty dictionary to make the csv further down

    for variant, genotype in genotypes.items():
        obs_ratio = compute_allele_frequency(genotype) #observed ratio
        if abs(obs_ratio - exp_ratio) > args.threshold:
            distorted.append(variant)
        if csv:
            obs_exp["obs"].append(obs_ratio) #ratio of observed variants
            obs_exp["exp"].append(exp_ratio) #ratio of expected variants 
    #Appending ratios to the dictionary

"""Emily suggestion: 
 If you change the default to None rather than False, and change 
 if csv to if type(csv) is str, then you can have the user specify the name of the resulting CSV 
 file instead of the default name I gave it.
"""
    if type(csv) is str:
        df_obs_exp = pd.DataFrame([obs_exp]) #naming the Pandas dataframe; edit so user can name 
        df_obs_exp.to_csv("obs_exp.csv", index = False) #writing the dataframe w/ ratios to a csv; edit so user can name

    write_results(distorted, output=args.output) #output is a CSV

"""Emily suggestion: 
 Write a function to plot the distorted regions, beginning with one 
 default plot. Later, it can be built out with different options 
 for different types of plots.
"""

def plot_detect(args,): #new function to plot the distorted regions
 #TO BE BUILT OUT 