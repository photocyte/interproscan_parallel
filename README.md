# interproscan_parallel
A Nextflow pipeline that wraps InterProScan , including making it parallelized &amp; doing some plotting

It is fairly bespoke to my personal troubleshooting and needs, and may have some external dependencies to other Nextflow scripts or resources that I haven't figure out how to put fully on github yet.

### Input
* A polypeptide FASTA file

### Running the workflow
See `run.sh` for how to run it. There are three workflows you might use via the `-entry` parameter for Nextflow:
* `download_data` that downloads the accessory data for the InterProScan Docker image
* The unnamed default workflow `` , that runs the InterProScan search, reformatting, and plotting. Best run on <10 polypeptides
* `do_simple_nonparallel_scan`, a simplified workflow that only does the InterProScan search and doesn't try to parallelize. Fine to run on 10,000+ polypeptides.


See `datasets/renaming.tsv` for how to structure this 2-column TSV. It allows taking the InterProScan member database unique identifiers, and maps them to a human readable name.

### Dependencies
* Nextflow (https://www.nextflow.io/)
* Apptainer (https://apptainer.org/)
* Miniconda (https://docs.conda.io/en/latest/) 

### Known issues
Despite the name, the parallelization is currently disabled though the vestigial code is mostly there and could be re-enabled if desired. The problem is, when InterProScan is run as a naive parallelization, the resulting GFF3 files are not mergable without breaking the specification (ID and/or Name attribute collisions). This renaming to avoid collisions and allow merging is partially implemented.  

