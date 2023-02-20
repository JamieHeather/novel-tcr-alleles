# Novel TCR alleles
## version 0.3.0 
### JH @ MGH, 2022

This repo aims to gather together published novel human TCR gene alleles as-yet not included in IMGT/GENE-DB, either inferred from full-variable domain spanning TCR-seq reads, observed directly in long-read genomic DNA, or from some similar technology that doesn't rely on assembling full regions from short reads. This may prove useful for TCRseq analysis applications, or for quick reference for efforts such as the [AIRR-C's Inferred Allele Review Committee](https://www.antibodysociety.org/the-airr-community/airr-subcomittees/inferred-allele-review-committee-iarc/), which aim to more rigorously update the field's germline knowledge.

### Input data

Currently this table collates alleles from:

* [Rodriguez *et al*., 2022, *BioRvix*](https://doi.org/10.1101/2022.05.24.493244)
* [Omer & Peres *et al*., 2022, *Genome Medicine*](https://doi.org/10.1186/s13073-021-01008-4)
* [Heather *et al*., 2022, *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkac190)
* [Lin *et al*., 2022, *Frontiers in Immunology*](https://doi.org/10.3389/fimmu.2022.922513)
* [Corcoran *et al*., 2023, *Immunity*](https://doi.org/10.1016/j.immuni.2023.01.026)

##### Notes on input data

* The differences in naming between the alleles found in both the Heather and Omer/Peres datasets is on account of me naming the ones I discovered using ungapped sequence positions, and them using gapped.
* The Corcoran *et al.* study includes five sets of monozygotic twins, which might affect consideration of the number of donors a given allele occurs in.  

### Compiling the data

The table can be generated using any combination of input files, provided they are supplied in the appropriate format and location. A single tab-separated file should be placed in the `input-data/` directory for each dataset to be included. Note that the name of these files will be used for column headers, so should be sensibly named. The fields that these files must contain are laid out in the `input-data/input-template.tsv` file. The fields, which should all have a value provided, are:

* Gene
  * The IMGT name of the gene the novel allele is from
* Novel-Allele-Name
  * The name or identifier for this allele in the originating dataset (if none provided just number them incrementally)
* Novel-Allele-Sequence
  * The nucleotide sequence of the novel allele, being either the V-REGION or J-REGION as defined by IMGT
* Number-Donors-With
  * The number of donors that were detected to contain this novel allele in this dataset (just label '1' if not known)
* Number-Donors-Searched
  * The total number of individuals that were screened in the assay where the allele was detected (label '1' if not known)

Running `compile-alleles.py` in the `scripts/` dir will read in all of the suitably provided input files, and concatenate them together, outputting a table where each row is a unique novel allele, with columns annotating the relevant details from the datasets in which they appear. This output file will appear in the root directory, named in the format `novel-TCR-alleles-YYYY-MM-DD-vX.X.X.tsv`, based on the date of running and the version of the repo, which should increment per changes to the input data. Any existing output files fitting that format in the root directory will be moved to the `archive/` directory.

Note that apart from selection of studies, no specific QC or additional validation has been applied to the sequences in question. The summary file does however contain columns of total numbers of datasets and donors each sequence was observed in.
