# Novel TCR alleles
## version 0.4.0 
### JH @ MGH, 2023

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

### Prerequisites

* `python3.10` or greater
  * `biopython`
  * `receptor_utils` (accessible in the user's PATH)

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

The table aims to:
* Combine the novel alleles presented in each input study into a single document, allowing a quick way to see how many of these studies each alleles occur in. Additionally, where possible, it aims to include information as to how many donors each allele was found in, and what identifiers they were referred to in each paper.
* Determine whether and when a novel allele may have been identified and published in an official IMGT release (potentially in a truncated form). This is done by comparing the sequences with the [historical entries of GENE-DB collated here](https://github.com/JamieHeather/genedb-releases).
* Assign a simple convenient name for the novel allele, in which it is named after the most closely related named IMGT allele, suffixed with information as to the polymorphisms which distinguish them. This is based on the [principles outlined here](https://wordpress.vdjbase.org/index.php/vdjbase_help/airr-seq-data-allele-names/) and achieved using William Lees' [`receptor_utils` package](https://williamdlees.github.io/receptor_utils/_build/html/index.html). This naming is calculated relative to the most recent IMGT/GENE-DB entry in the `genedb-releases` submodule, and therefore may be subject to change.
  * Note that this process does not work for a small number of alleles for which there isn't a suitable reference in the IMGT release used, prompting a `No reference genes to compare` error from `receptor_utils`. This mostly seems to occur for out-of-frame pseudogenes, which gapped sequences are not available for.
* The `receptor_utils` package also provides other useful information, such as IMGT gapped sequences and warnings when novel alleles are missing conserved residues; this information is also retained in the output table.

Running `compile-alleles.py` in the `scripts/` dir will read in all of the suitably provided input files, and concatenate them together, outputting a table where each row is a unique novel allele, with columns annotating the relevant details from the datasets in which they appear. This output file will appear in the root directory, named in the format `novel-TCR-alleles_YYYY-MM-DD_vX.X.X_GENEDB-XXXXX-X.tsv`, based on the date of running, the version of the repo (which should increment per changes to the input data or table generation code), and the [IMGT/GENE-DB release](https://www.imgt.org/download/GENE-DB/) used as a comparison. Any existing output files fitting that format in the root directory will be moved to the `archive/` directory.

Note that apart from selection of studies, no specific QC or additional validation has been applied to the sequences in question. The summary file does however contain columns of total numbers of datasets and donors each sequence was observed in.
