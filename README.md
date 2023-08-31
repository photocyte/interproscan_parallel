# interproscan_parallel
A Nextflow pipeline that wraps InterProScan , including making it parallelized &amp; doing some plotting

It is fairly bespoke to my personal troubleshooting and needs, and has some external dependencies to other Nextflow scripts that I haven't figure out how to put fully on github yet.

But in principle, should be able to run it with this:
```
nextflow run https://github.com/photocyte/interproscan_parallel/ -resume --fasta polypeptides_to_annotate.pep.fa
```
