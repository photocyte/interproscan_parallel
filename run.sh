#NF_SCR="https://github.com/photocyte/interproscan_parallel"
NF_SCR="$(pwd)/interproscan_parallel.nf"

## Downloading the accessory data for Docker/Singularity/Apptainer execution of InterProScan
## This ends up in the dependencies/download_interproscan_docker_data/ folder.
## You might have to modify the download_interproscan_docker_data process, and the nextflow.config singularity.runOptions to get it to work.
##mkdir -p dependencies/download_interproscan_docker_data/interproscan-5.73-104.0/

### Setting up the various TMPDIRs . Use the /dev/shm ramdisk if performance is a concern 
TOTAL_TMPDIR="/dev/shm/$USER/container_executor_cache" 
#TOTAL_TMPDIR="$(pwd)/container_executor_cache" 
mkdir -p "${TOTAL_TMPDIR}"
#export APPTAINER_CACHEDIR="${TOTAL_TMPDIR}"
#export APPTAINER_TMPDIR="${TOTAL_TMPDIR}"
export SINGULARITY_CACHEDIR="${TOTAL_TMPDIR}"
export SINGULARITY_TMPDIR="${TOTAL_TMPDIR}"
export TMPDIR="${TOTAL_TMPDIR}"

#nextflow run ${NF_SCR} -c conf/download.config -latest -resume -entry download_data ## Used to download the interproscan docker/singularity dependency data (HMMs, etc.). It is about 50 GB in total.

## An example of how to symlink the accessory data, if interproscan_parallel was run as a subdirectory of the parent directory:
#mkdir -p dependencies/download_interproscan_docker_data/interproscan-5.73-104.0/
#cd dependencies/download_interproscan_docker_data/interproscan-5.73-104.0/
#ln -s ../../../../interproscan_parallel/dependencies/download_interproscan_docker_data/interproscan-5.73-104.0/data/ .
#cd ../../../

## Running the InterProScan annotation and graphical plotting for the peptide FASTA of the 3 PKZILLAs:
nextflow run ${NF_SCR} -r main -latest -resume -profile ramdisk --fasta PKZILLA-1-on-B.loci.nt.rename.consensus.fa.pep.fa

