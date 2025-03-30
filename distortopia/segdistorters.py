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
            obs_exp["obs"].append(obs_ratio)
            obs_exp["exp"].append(exp_ratio) #Appending ratios to the dictionary

    if csv:
        df_obs_exp = pd.DataFrame([obs_exp])
        df_obs_exp.to_csv("obs_exp.csv", index = False) #Writing the ratios to a csv

    write_results(distorted, output=args.output)
