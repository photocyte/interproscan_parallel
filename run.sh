#ln -s ../../20230726_12B1_v1.1/gene_models/12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa.gz .
#seqkit grep -n -r -p "PKZILLA" 12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa.gz > PKZILLA_12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa

NF_SCR="https://github.com/photocyte/interproscan_parallel"
#NF_SCR=/home/tfallon/source/nextflow/interproscan_parallel/interproscan_parallel.nf

## Downloading the accessory data for Docker/Singularity execution of InterProScan
## This ends up in the dependencies/download_interproscan_docker_data/ folder.
## You might have to modify the download_interproscan_docker_data process, and the nextflow.config singularity.runOptions to get it to work.
##mkdir dependencies/download_interproscan_docker_data/interproscan-5.64-96.0/data
#nextflow run ${NF_SCR} -latest -resume -entry download_data ## Used to download the interproscan docker data thingy

## Running the InterProScan annotation and graphical plotting for the peptide FASTA of the 3 PKZILLAs:
nextflow run ${NF_SCR} -r main -latest -resume --fasta PKZILLA-1-on-B.loci.nt.rename.consensus.fa.pep.fa --renaming ./dataset/renaming.tsv

