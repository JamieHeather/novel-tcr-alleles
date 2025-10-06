# Novel TCR alleles
## version 0.7.0 
### JH @ MGH, 2023

This repo aims to gather together published novel human TCR gene alleles as-yet not included in IMGT/GENE-DB, either inferred from full-variable domain spanning TCR-seq reads, observed directly in long-read genomic DNA, or from some similar technology that doesn't rely on assembling full regions from short reads. This may prove useful for TCRseq analysis applications, or for quick reference for efforts such as the [AIRR-C's Inferred Allele Review Committee](https://www.antibodysociety.org/the-airr-community/airr-subcomittees/inferred-allele-review-committee-iarc/), which aim to more rigorously update the field's germline knowledge.

### Input data

Currently this table collates alleles from:

* [Rodriguez *et al*., 2022, *Cell Genomics*](https://doi.org/10.1016/j.xgen.2022.100228)
* [Omer & Peres *et al*., 2022, *Genome Medicine*](https://doi.org/10.1186/s13073-021-01008-4)
* [Heather *et al*., 2022, *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkac190)
* [Lin *et al*., 2022, *Frontiers in Immunology*](https://doi.org/10.3389/fimmu.2022.922513)
* [Corcoran *et al*., 2023, *Immunity*](https://doi.org/10.1016/j.immuni.2023.01.026)
* [Mikelov *et al*., 2024, *Genome Research*](https://doi.org/10.1101/gr.278775.123)
* [Mantena *et al*., 2025, *bioRxiv*](https://doi.org/10.1101/2025.08.20.671277 )

##### Notes on input data

* The differences in naming between the alleles found in both the Heather and Omer/Peres datasets is on account of me naming the ones I discovered using ungapped sequence positions, and them using gapped.
* The Corcoran *et al.* study includes five sets of monozygotic twins, which might affect consideration of the number of donors a given allele occurs in.  
* The Mikelov *et al.* data are described in the paper linked above, but actually come from the associated [VDJ.online 'Gene Library'](https://vdj.online/) resource (accessed on 2023-10-16).
  * Note that the donor count information for this dataset is actually technically a haplotype count. 
* The Mantena *et al.* data described in the above linked pre-print actually come from the [associated GitHub repo (accessed on 2025-10-05)](https://github.com/SreekarMantena/tcrdiversity). 
  * Note that these data did not come with donor count numbers, therefore each allele was arbitrarily said to be found in 1 donor

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

Running `compile-alleles.py` in the `scripts/` dir will read in and integrate all of the suitably provided input files,outputting a table where each row is a unique novel allele, with columns annotating the relevant details from the datasets in which they appear. This output file will appear in the root directory, named in the format `novel-TCR-alleles_YYYY-MM-DD_vX.X.X_GENEDB-XXXXX-X.tsv`, based on the date of running, the version of the repo (which should increment per changes to the input data or table generation code), and the [IMGT/GENE-DB release](https://www.imgt.org/download/GENE-DB/) used as a comparison. Any existing output files fitting that format in the root directory will be moved to the `archive/` directory. 

The output table aims to:
* Combine the novel alleles presented in each input study into a single document, allowing a quick way to see how many of these studies each alleles occur in. Additionally, where possible, it aims to include information as to how many donors each allele was found in, and what identifiers they were referred to in each paper.
* Determine whether and when a novel allele may have been ...:
  * Published by IMGT in an official GENE/DB release (potentially in a truncated form). 
    * This is done by comparing the sequences with the [historical entries of GENE-DB collated here](https://github.com/JamieHeather/genedb-releases).
  * Affirmed by IARC.
    * This is achieved by checking [the appropriate page of the Open Germline Receptor Database (OGRDB)](https://ogrdb.airr-community.org/sequences/Human_TCR). 
* Assign a simple convenient name for the novel allele, in which it is named after the most closely related named IMGT allele, suffixed with information as to the polymorphisms which distinguish them. This is based on the [principles outlined here](https://wordpress.vdjbase.org/index.php/vdjbase_help/airr-seq-data-allele-names/) and achieved using William Lees' [`receptor_utils` package](https://williamdlees.github.io/receptor_utils/_build/html/index.html). This naming is calculated relative to the most recent IMGT/GENE-DB entry in the `genedb-releases` submodule, and therefore may be subject to change.
  * Note that this process does not work for a small number of alleles for which there isn't a suitable reference in the IMGT release used, prompting a `No reference genes to compare` error from `receptor_utils`. This mostly seems to occur for out-of-frame pseudogenes, which gapped sequences are not available for.
  * The `receptor_utils` package also provides other useful information, such as IMGT gapped sequences and warnings when novel alleles are missing conserved residues; this information is also retained in the output table.
* Find deposited matches and potential orthogonal validation of novel sequences.
  * The complete gene sequence is then BLASTed against NCBI's nt database using the BioPython 'NCBIWWW' function, to check for existing deposited data that may help validate or confirm. 
  * This retains only 100% matches to human sequences.
  * Note that this process can take a long time, and thus users may wish to comment it out if they're running it themselves and don't need this data.
  * Also note that some of the validations of these alleles by the depositing groups have been uploaded to GenBank, so some of the accessions reported relate to evidence from the same original paper.

Note that apart from selection of studies, no specific QC or additional validation has been applied to the sequences in question. The summary file does however contain columns of total numbers of datasets and donors each sequence was observed in.

#### Related citation

This resource can be cited via [this 2025 *Immunoinformatics* manuscript I wrote with some AIRR-C colleagues, "*The gremlin in the works: why T cell receptor researchers need to pay more attention to germline reference sequences*"](10.1016/j.immuno.2025.100058).
