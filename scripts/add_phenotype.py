#! /usr/bin/env python

import requests, sys, os
import csv
import argparse
import gzip, shutil

parser = argparse.ArgumentParser(description="Add phenotype to outrider output", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("file", help="provide outrider input file to add phenotype")
parser.add_argument("-t", "--type", type=int, help="Result file input type. 0 for outrider, 1 for featurecounts, 2 for fraser")
parser.add_argument("-p", "--pval", type=float, help="Max p-value to search phenotype for")
parser.add_argument("-z", "--zscore", type=float, help="Min z-score to search phenotype for")
parser.add_argument("-m", "--meancorr", type=int, help="Min mean corrected to search phenotype for")

def get_phenotype(gene):
    '''Returns ensemble phenotype'''
    server = "https://rest.ensembl.org"
    ext = "/phenotype/gene/homo_sapiens/"
    output_pheno=""
    output_omim=""
    # try:
    r = requests.get(server+ext+gene, headers={ "Content-Type" : "application/json"})
    # except:
    if not r.ok:
        output_pheno = "error in request"#r.raise_for_status()
        # sys.exit()
    else:
        decoded = r.json()
        for i in decoded:
            if i['source']=='MIM morbid':
                output_omim += i['source']+':'+i['description']+";" 
            output_pheno += (i['source']+':'+i['description']+";")
        if output_omim: output_pheno = output_omim
        elif not output_pheno: output_pheno = "nothing found"
    return output_pheno.replace(",","-").replace('\"',"")

def return_pheno_row(type_res, row, output_pheno):
    '''returns new row including added phenotype
    for each type of input result file'''
    if type_res==1: 
        row.insert(-2,output_pheno)
    elif type_res==0:
        row.append(output_pheno)
    elif type_res==2:
        row.insert(len(row)+1,output_pheno)
    return(row)

def add_pheno_to_row(type_res,i,row,header_row,pval_bound,zscore_bound,meancorr_bound):
    '''get phenotype
    add as new column to row
    and return as newrow'''
    output_pheno="outside boundaries"
    if i==0:
        output_pheno="phenotype"
    #for feature count input
    elif type_res==1: 
        output_pheno = get_phenotype(row[header_row.index('gene_name')])
    #for outrider result input
    elif type_res==0 and abs(float(row[header_row.index('zScore')])) > zscore_bound and float(row[header_row.index('pValue')]) < pval_bound and float(row[header_row.index('meanCorrected')]) >= meancorr_bound:
        output_pheno = get_phenotype(row[header_row.index('gene')])
    #for fraser result input
    elif type_res==2  and float(row[header_row.index('pValue')+1]) < pval_bound:
        # and float(row[header_row.index('meanCorrected')]) >= meancorr_bound
        gene = row[header_row.index('hgncSymbol')+1]
        if ";" in gene: 
            output_pheno = ""
            for g in gene.split(";"):
                output_pheno += get_phenotype(g) + ";" 
        else: output_pheno = get_phenotype(gene)
    
    return(return_pheno_row(type_res,row,output_pheno))
    
def add_phenotype_to_file(file_in, type_res=0, pval_bound=0.01, zscore_bound=3, meancorr_bound=0):
    '''Add phenotype as a column to imported file
    for rows that match boundaries on pvalue, zscore and/or meancorr value
    and write in new file'''
    file_out = os.path.splitext(file_in)[0]+"_pheno_added"+os.path.splitext(file_in)[1]
    #file_out = file_in.split(".")[0]+"_pheno_added."+ext
    delim = ","
    if type_res != 0: delim = "\t"
    with open(file_out, 'w') as newcsvfile:
        writer = csv.writer(newcsvfile, delimiter=delim)
        with open(file_in, newline='') as csvfile:
            reader = csvfile.readlines()
            header_row = reader[0][0:-1].split(delim)
            for i,row in enumerate(reader):
                newrow = add_pheno_to_row(type_res,i,row[0:-1].split(delim),header_row,pval_bound,zscore_bound,meancorr_bound)
                writer.writerow(newrow)
    if type_res == 0:
        with open(file_out, 'rb') as f_out:
            with gzip.open(file_out+".gz", 'wb') as f_outgz:
                shutil.copyfileobj(f_out, f_outgz)
    print("done")

def main():
    '''Set and get arguments to add phenotype to results within corresponding boundaries'''
    args = parser.parse_args()
    file = args.file
    type_res = 0 
    pval = 0.01
    zscore = 2.5
    meancorr = 0
    if args.type is not None: type_res = args.type
    if args.pval is not None: pval = args.pval
    if args.zscore is not None: zscore = args.zscore
    if args.meancorr is not None: meancorr = args.meancorr
    
    add_phenotype_to_file(file, type_res, pval, zscore, meancorr)

main()
