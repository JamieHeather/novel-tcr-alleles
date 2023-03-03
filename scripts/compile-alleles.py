# -*- coding: utf-8 -*-

"""
compile-alleles

Go through all the TSV files in input-data/, and compile into one large summary file
"""

import os
import collections as coll
import pandas as pd
import datetime
import subprocess

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.3.0'
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


def read_fa(ff):
    """
    :param ff: opened fasta file
    read_fa(file):Heng Li's Python implementation of his readfq function (tweaked to only bother with fasta)
    https://github.com/lh3/readfq/blob/master/readfq.py
    """

    last = None                                 # this is a buffer keeping the last unprocessed line
    while True:                                 # mimic closure
        if not last:                            # the first record or a record following a fastq
            for l in ff:                        # search for the start of the next record
                if l[0] in '>':                 # fasta header line
                    last = l[:-1]               # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for l in ff:                            # read the sequence
            if l[0] in '>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':          # this is a fasta record
            yield name, ''.join(seqs), None     # yield a fasta record
            if not last:
                break
        else:
            raise IOError("Input file does not appear to be a FASTA file - please check and try again")


def flatten(list_of_lists):
    """
    :param list_of_lists: List of lists, as may be generated at certain parts of the analysis
    :return: Said list of lists flattened into a single one-level deep list
    """
    return [thing for sublist in list_of_lists for thing in sublist]


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

# If that all completed successfully, move any existing compilations from root folder to the archive
out_prefix = 'novel-TCR-alleles-'
older_versions = [x for x in os.listdir('../') if x.startswith(out_prefix) and
                  (x.endswith('.tsv') or x.endswith('.tsv.gz'))]

currently_archived = os.listdir('../archive/')
for ov in older_versions:
    if ov not in currently_archived:
        os.replace('../' + ov, '../archive/' + ov)
    else:
        os.replace('../' + ov, '../archive/' + ov.replace('.tsv', '-name-clash.tsv'))

# Then stick it all together into a table before further processing
out_dat = list_to_df(out_dat, out_headers, False)

# Need to determine whether the alleles are featured in IMGT/GENE-DB, and if so when
# Go through the genedb-releases scrapes, and compile allele-specific release data
release_path = '../genedb-releases/releases/'
release_dirs = [x for x in os.listdir(release_path) if os.path.isdir(release_path + x) and x.startswith('20')]
release_dirs.sort()
genedb = coll.defaultdict(list)
partial = coll.defaultdict()
for rdir in release_dirs:
    release = rdir.split('_')[-1]
    genedb_file = [x for x in os.listdir(release_path + rdir) if 'fasta-nt-WithoutGaps-F+ORF+allP' in x]
    if len(genedb_file) != 1:
        raise IOError('Unexpected number of file matches in release folder ' + rdir)
    else:
        with open(release_path + rdir + '/' + genedb_file[0], 'r') as in_file:
            for header, seq, qualnull in read_fa(in_file):
                # Take human TCR V/J genes
                if 'Homo sapiens' not in header:
                    continue
                bits = header.split('|')
                if bits[1][:2] != 'TR' or bits[1][3] not in ['V', 'J']:
                    continue
                seq = seq.upper()
                genedb[seq].append(bits[1] + '|' + release)
                if bits[13] != ' ':
                    partial[bits[1]] = bits[13].strip()

# Then go back through our novel alleles, and see if any appear in a previous release
# Simultaneously also use William's name_allele utility to generate a standardised name, using most recent IMGT release

# First generate the necessary reference files
gapped_file = release_path + rdir + '/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'
for locus in ['A', 'B', 'G', 'D']:
    cmd = 'extract_refs -L TR' + locus + ' ' + gapped_file + ' "Homo sapiens"'
    subprocess.call(cmd, shell=True)

# Then iterate across the dataframe
imgt_release = []
imgt_id = []
std_id = []
gapped = []
notes = []
file_suffix = {'V': '_gapped.fasta ', 'J': '.fasta '}

# Check whether each 'novel' allele is actually an existing covered IMGT gene
for na in out_dat.index:
    # Check for previous appearances
    seq = out_dat.loc[na]['Ungapped-Sequence']
    tmp_imgt_id, tmp_imgt_release, note, ids = [''] * 4

    # First check for exact matches ...
    if seq in genedb:
        ids = list(set([x.split('|')[0] for x in genedb[seq]]))
        releases = list(set([x.split('|')[1] for x in genedb[seq]]))
        note += 'Exact match of '

    # ... Then check whether the novel allele is a substring of any IMGT allele ...
    elif len([x for x in genedb if seq in x]) > 0:
        ids = list(set([x.split('|')[0] for x in flatten([genedb[x] for x in genedb if seq in x])]))
        releases = list(set([x.split('|')[1] for x in flatten([genedb[x] for x in genedb if seq in x])]))
        note += 'Shorter version of '

    # ... Or whether any IMGT allele is a substring of it
    elif len([x for x in genedb if x in seq]) > 0:
        ids = list(set([x.split('|')[0] for x in flatten([genedb[x] for x in genedb if x in seq])]))
        releases = list(set([x.split('|')[1] for x in flatten([genedb[x] for x in genedb if x in seq])]))
        note += 'Longer version of '

    # Account for potential multiple matches
    if ids:
        ids.sort()
        tmp_imgt_id = '|'.join(ids)
        note += tmp_imgt_id
        if tmp_imgt_id in partial:
            note += ' (' + partial[tmp_imgt_id] + '). '
        else:
            note += '. '

        releases.sort()
        tmp_imgt_release = releases[0]

    imgt_release.append(tmp_imgt_release)
    imgt_id.append(tmp_imgt_id)

    # Then generate standardised naming (off this release)
    gene = out_dat.loc[na]['Gene']
    cmd = 'name_allele -g ' + gene + ' Homo_sapiens_' + gene[:4] + file_suffix[gene[3]] + seq

    naming = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL).decode('UTF-8').rstrip()
    if '\n' in naming:
        naming_bits = naming.split('\n')
        if len(naming_bits) in [3, 5]:
            std_id.append(naming_bits[0].upper())
            gapped.append(naming_bits[2])
            if len(naming_bits) == 5:
                note += naming_bits[4] + '. '
        else:
            raise IOError("Unknown format")
    else:
        gapped.append('')
        if naming == 'No reference genes to compare':
            std_id.append('-')
        else:
            std_id.append(naming.upper())

    notes.append(note)

# Add those new columns back in, and write the updated table out
out_dat.insert(1, 'IMGT-ID', imgt_id)
out_dat.insert(1, 'IMGT-Appearance', imgt_release)
out_dat.insert(1, 'Standard-ID', std_id)
out_dat.insert(5, 'Gapped-Sequence', gapped)
out_dat.insert(4, 'Notes', notes)

out_nam = out_prefix + '_'.join([today(), 'v' + run_version, 'GENEDB-' + release + '.tsv'])
out_dat.to_csv('../' + out_nam, sep='\t', index=False)

# ... and tidy up the references
subprocess.call('rm Homo*fasta', shell=True)
