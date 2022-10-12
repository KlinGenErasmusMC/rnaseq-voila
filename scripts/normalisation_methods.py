#! /usr/bin/env python3


"""
Read count normalisation with TPM; TMM;

At the moment per exon
"""



import os
import argparse
import sys
import pandas as pd
import numpy as np
import logging
import subprocess


R_TMM = os.path.join( os.path.dirname(os.path.realpath(__file__)), "tmm.R")
 
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, help="input file with read counts for all samples", required=True)
    parser.add_argument("-t", "--tpm-output", type=str, help="output file to write TPM to")
    parser.add_argument("-m", "--tmm-output", type=str, help="output file to write TMM to, requires R in PATH with libraries argparse and edgeR installed")
    parser.add_argument("-l", "--logfile", type=str, default="{}.log".format(os.path.basename(__file__)), help="Name of the log file [{}]".format(os.path.basename(__file__)))
    parser.add_argument("-a", "--aggregate", action='store_true', help="aggregate per gene")
    
    args = parser.parse_args()

    logging.basicConfig(filename=args.logfile, level=logging.DEBUG)
    logging.info("Running {} with the followig options: {}".format(__file__,args))
    
    # read counts
    df = pd.read_csv(args.input, sep="\t")

    samples = df.columns[4:]

    # calculate length of items
    df['length'] = df['end'] - df['start']
    df[df['length'] == 0] = 1  # fix for in case length is 0
    columns = list(df)

    # aggregate
    if args.aggregate:
        df = df.groupby(["chr", "gene"]).agg("sum").reset_index().reindex(columns=columns)
           
    # get flagstats + correction factor)

    if args.tpm_output:
        TPM = df.copy()
        for sample in samples:
           A = TPM[sample] * 10**3 / TPM['length']
           TPM[sample] = A * (1/A.sum()) * 10**6 if A.sum() > 0 else 0
        TPM.drop(df.columns[len(df.columns)-1], axis=1, inplace=True)
        TPM.drop(columns=["start", "end"])
        TPM.to_csv(args.tpm_output, sep="\t", index=False)

    if args.tmm_output:
        if args.aggregate:
            agg = df.drop(df.columns[len(df.columns)-1], axis=1, inplace=False).to_csv(args.input + ".agg", sep="\t", index=False)
            f_in = args.input + ".agg"
        else:
            f_in = args.input
        command = [R_TMM, "--count-table", f_in, "--outfile", args.tmm_output]
        subprocess.call(command, stdout=subprocess.DEVNULL)
