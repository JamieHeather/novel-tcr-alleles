# -*- coding: utf-8 -*-

"""
compile-alleles

Go through all the TSV files in input-data/, and compile into one large summary file
"""

import os
import collections as coll
import pandas as pd
import datetime

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.0'
__author__ = 'Jamie Heather'


def nest():
    """
    :return: a list defaultdict, for nesting multiple layers of defaultdicts
    """
    return coll.defaultdict(list)


def list_to_df(input_list, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param input_list: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(input_list)
    df = df.rename(index=str, columns=dict(zip(range(len(headers)), headers)))
    df = df.sort_values(by=[headers[0]])
    if rename is True:
        df = df.set_index(headers[0], drop=True)
    return df


def today():
    """
    :return: Today's day, in ISO format
    """
    return datetime.datetime.today().date().isoformat()


# Detect data files
data_dir = '../input-data/'
input_files = [x for x in os.listdir(data_dir) if x.endswith('.tsv') and 'template' not in x]
input_files.sort()

# Figure out version number
with open('../README.md', 'r') as in_file:
    for line in in_file:
        if line.startswith('## version'):
            run_version = line.rstrip().split(' ')[-1]
            break

if 'run_version' not in locals():
    raise IOError("Cannot determine version number, check README formatting.")

# Go through and read each in, also determining unique novel alleles by sequence
uniq_seqs = coll.defaultdict(list)
dfs = coll.defaultdict()
dicts = coll.defaultdict()
for fl in input_files:
    nam = fl[:-4]
    dfs[nam] = pd.read_csv(data_dir + fl, sep='\t')
    dicts[nam] = coll.defaultdict()
    for row in dfs[nam].index:
        row_dat = dfs[nam].loc[row]
        seq = row_dat['Novel-Allele-Sequence'].upper()
        uniq_seqs[seq].append(row_dat['Gene'])
        dicts[nam][seq] = coll.defaultdict(dict, {
            'Gene': row_dat['Gene'],
            'Novel-Allele-Name': row_dat['Novel-Allele-Name'],
            'Number-Donors-With': row_dat['Number-Donors-With'],
            'Number-Donors-Searched': row_dat['Number-Donors-Searched']})


# Check the genes associated with each sequence are all the same
uniq_check = [x for x in uniq_seqs if len(list(set(uniq_seqs[x]))) > 1]
if uniq_check:
    for clash in uniq_check:
        print(clash, uniq_seqs[clash])
    raise IOError("Allele sequence(s) found with >1 gene associated - review above sequences/genes and try again.")

# Then go through these data and compile into a summary dataframe
sequences = list(uniq_seqs.keys())
datasets = list(dicts.keys())
datasets.sort()

# Set up the columns required, including the dataset-specific ones
out_headers = ['Gene', 'Ungapped-Sequence']
for dataset in datasets:
    out_headers.append(dataset + '-Name')
    out_headers.append(dataset + '-Donor-Count')
    out_headers.append(dataset + '-Donor-Total')

out_dat = []
for novel_seq in sequences:
    row_list = [uniq_seqs[novel_seq][0], novel_seq]
    number_datasets = 0
    number_donors = 0
    number_tested = 0
    for dataset in datasets:
        if novel_seq in dicts[dataset]:
            tmp = dicts[dataset][novel_seq]
            allele_nam = tmp['Novel-Allele-Name']
            donors_with = tmp['Number-Donors-With']
            donors_tested = tmp['Number-Donors-Searched']
            row_list += [allele_nam, donors_with, donors_tested]
            number_datasets += 1
            number_donors += donors_with
            number_tested += donors_tested
        else:
            row_list += [''] * 3
    out_headers += ['Number-Datasets-In', 'Number-Donors-In', 'Number-Donors-Tested']
    row_list += [number_datasets, number_donors, number_tested]
    out_dat.append(row_list)

# Then stick it all together and write it out
out_dat = list_to_df(out_dat, out_headers, False)
out_nam = 'novel-TCR-alleles-' + today() + '-v' + run_version + '.tsv'

out_dat.to_csv('../' + out_nam, sep='\t', index=False)
