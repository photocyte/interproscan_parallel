// Also, see here https://lipm-gitlab.toulouse.inra.fr/LIPM-BIOINFO/nextflow-interproscan5/blob/master/interproscan.nf
nextflow.enable.dsl=2

//This isn't the most elegant way to split up the modules.
//It's currently not possible in Nextflow (https://github.com/nextflow-io/nextflow/discussions/4254), and may not ever be implemented (https://github.com/nextflow-io/nextflow/issues/4112)
include { generate_svg_colors ; svg_2_pdf ; svg_utils_merge ; svg_resize_page ; split_gff_by_seqid ; DNA_features_viewer ; pdf_2_PDF_A_1B } from './modules/dna_features_viewer/dna_features_viewer.nf'

//Seeing as the unnamed main workflow, and do_simple_nonparallel_scan, currently do not call these processes, commenting them out. 
//In short, this was a fairly complicated way to to target selection, within the script. The alternative (and what I'm doing now) is to use some simple grepping outside the script
//include { detab ; debar ; decolon ; dummy_publish_path } from '/home/tfallon/source/nextflow/nf_utility/nf_utility.nf'
//include { download_pfam_A ; prepare_pfam_A ; pfam_run ; hmmpress_generic ; hmm_subset } from '/home/tfallon/source/nextflow/pfam_nf/pfam.nf'

process download_ipr2go {
storeDir "results/${task.process}"
output:
 path "interpro2go.txt"
shell:
'''
wget --no-check-certificate http://www.geneontology.org/external2go/interpro2go -O interpro2go.txt
'''
}

process download_interproscan_docker_data {
storeDir "dependencies/download_interproscan_docker_data"
output:
 path "interproscan-*/data"
shell:
'''
##Be sure to check the 'nextflow.config' for how to mount this 'data' directory into the ultimate interproscan_run process

VERSION='5.73-104.0' ##formerly 5.64-96.0
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/alt/interproscan-data-${VERSION}.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/alt/interproscan-data-${VERSION}.tar.gz.md5
md5sum -c interproscan-data-${VERSION}.tar.gz.md5
tar -pxzf interproscan-data-${VERSION}.tar.gz

## Files & executables are supplied by the Docker image.
ln -s /opt/interproscan/interproscan.properties interproscan-*/
ln -s /opt/interproscan/setup.py interproscan-*/
ln -s /opt/interproscan/bin interproscan-*/

cd interproscan-*/
python3 setup.py -f interproscan.properties ## This does hmmpress, and some other things. Not mentioned in the InterProScan Docker image documentation, but I think it's important.
cd ../

##Clean up symlinks & tar.gz file
rm -f interproscan-*/interproscan.properties
rm -f interproscan-*/setup.py
rm -f interproscan-*/bin
rm -f interproscan-data-${VERSION}.tar.gz
'''
}

process interproscan_run {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    cpus 8
    //scratch true
    //scratch 'ram-disk' // Can't use this, plus nextflow.config singularity.runOptions -B mounting at the same time.
    //stageInMode 'copy'
    errorStrategy 'finish'
    maxRetries 3
    executor "local"
    queue "home"
    time '24h'
    maxForks 10
    cache "deep"
    debug true //Print stdout to the terminal
    input:
     path inputFasta // Might consider linking in the nextflow.config singularity.runOptions -B 'data' directory via a Channel instead.
    output:
     path "*.gff3", emit: gffs
     path "*.json", emit: jsons
     path "*.tsv", emit: tsvs 
     path "*.tsv.sites", emit: sites, optional:true
     path "*.xml", emit: xmls
    shell:
      '''
      ##Old way, using the manually installed interproscan 
      ##/home/tfallon/software/interproscan-5.63-95.0/interproscan.sh -appl ${APPLS} --disable-precalc --iprlookup --pathways --goterms --enable-tsv-residue-annot --cpu !{task.cpus} -i !{inputFasta} --tempdir ./

      ##Removed MobiDBLite, as it had been crashing
      ##CDD was also crashing with a rpsbproc dependency issue
      ##APPLS="COILS,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM"
      ##Will try running with all tools, since it is now the Docker container. Readd this parameter if wanting to disable again: '-appl ${APPLS}' 

      ## New way, using a Singularity/Apptainer execution off the official interproscan Docker container
      /opt/interproscan/interproscan.sh --disable-precalc --iprlookup --pathways --goterms --enable-tsv-residue-annot --cpu !{task.cpus} -i !{inputFasta} --tempdir ./
      '''
}

process hmm_subset {
input:
 val target_hmm
 path hmmfile
output:
 path "${target_hmm}.hmm"
tag "${target_hmm}"
shell:
'''
cat !{hmmfile} | python -c "import re;import sys;data = sys.stdin.read(); [ print(m) for m in re.findall('HMMER3.{1,50}?ACC[\s]{1,6}!{target_hmm}.+?^//',data,flags=re.DOTALL|re.MULTILINE)]" | tee !{target_hmm}.hmm
'''
}

process fasta_remove_asterisk {
cpus 1
input:
 path fasta
output:
 path "results/${fasta}"
shell:
'''
##Note the additional backslash for Nextflow escaping. Replace internal '*'s with Xs.
mkdir results
seqkit replace -s -p "\\*$" -r "" !{fasta} | seqkit replace -s -p "\\*" -r "X" > results/!{fasta}
'''
}

process gff_strip_fasta {
cpus 1
input:
 path gff
output:
 path "results/ipr.${gff}", emit: gffs, optional: true
 path "results/ipr.${gff}.fasta", emit: fastas, optional: true
shell:
'''
mkdir results

GFF_LEN=`wc -l !{gff} | cut -f 1 -d " "`
if [ $GFF_LEN -eq 4 ]
then
##Empty GFF file, so skip it.
exit 0
fi

##Note: Below ID renaming only works for if it is 1 FASTA target per GFF
SEQID=`cat !{gff} | grep -v "##" | head -n 1 | cut -f 1` 
##TARGETS="##|Gene3D|polypeptide|SUPERFAMILY"

##In the past, this process also did filtering. with | grep -P ${TARGETS} > results/ipr.!{gff}

head -n +`grep -n "##FASTA" !{gff} | cut -d : -f 1` !{gff} | grep -v "##FASTA" > results/ipr.!{gff}
tail -n +`grep -n "##FASTA" !{gff} | cut -d : -f 1` !{gff} | grep -v "##FASTA" > results/ipr.!{gff}.fasta
'''
}

process gff_append_name_to_seqid_ipr {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cpus 1
input:
 path gff
output:
 path "results/${gff}"
shell:
'''
gff_append_name_to_seqid_ipr.py !{gff}
'''
}


process hmmersearch2gff {
input:
 path inputTbl
output:
 path "${inputTbl}.gff"
script:
"""
#!/usr/bin/env python
## Note the additional backslash for Nextflow escaping
import re
import os
r_handle = open("${inputTbl}","r")
w_handle = open("${inputTbl}.gff","w")
domainIndex = dict()
for l in r_handle.readlines():
    if l[0] == "#":
        continue
    sl = re.split('\s+', l)
    print(sl)
    pfam_domain_name = sl[3]
    pfam_domain_accession = sl[4]
    if pfam_domain_accession == "-":
         pfam_domain_accession = sl[3]
    seqid, source, type, start, end, score, strand, phase, attributes = sl[0], "PFAM", "hmmsearch", sl[19], sl[20], sl[11], "+", ".", ["ID="+pfam_domain_accession,"Name="+pfam_domain_name]
    if seqid not in domainIndex.keys():
         domainIndex[seqid] = {pfam_domain_name:1}
    else:
         if pfam_domain_name in domainIndex[seqid].keys():
             domainIndex[seqid][pfam_domain_name] += 1
         else:
              domainIndex[seqid][pfam_domain_name] = 1
    attributes[0] = "ID="+seqid+"_"+attributes[0][3:]
    attributes[0] += "-dom"+str(domainIndex[seqid][pfam_domain_name]) ## Add an index of domain number to the ID
    w_handle.write('\t'.join([seqid, source, type, start, end, score, strand, phase, ";".join(attributes)])+os.linesep)
r_handle.close()
w_handle.close()    
"""
}

process gff_score_filter {
input:
 path inputGff
 val scoreToFilterBy //e.g. 1E-3
output:
 path "filtered/${inputGff}"
conda 'bioawk genometools-genometools'
shell:
'''
mkdir filtered
cat !{inputGff} | bioawk -c gff '$score < !{scoreToFilterBy}{ print $0 }' | gt gff3 -tidy -sort -retainids > filtered/!{inputGff}
'''
}


process mgkit_hmmer2gff {
//conda "mgkit" //Will use docker instead... Conda too slow
//Turns out the Singularity container of the docker fails with an error.
input:
 tuple path(inputFasta),path(hmmsearch_tbl)
output:
 path("${hmmsearch_tbl}.gff")
shell:
'''
python ../../../pfam2gff.py -i !{hmmsearch_tbl} -e 1e-3 > !{hmmsearch_tbl}.gff
##perl ../../../extract_genes_with_domain.pl !{hmmsearch_tbl} !{hmmsearch_tbl}.gff
##hmmer2gff -o !{hmmsearch_tbl}.gff !{inputFasta} !{hmmsearch_tbl}
'''
}

process gff_nested_filter {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
conda 'bedops biopython'
input:
 path inputGff
output:
 path "filtered/${inputGff}", emit: remaining_gffs
 path "nested.gff", emit: removed_gffs
shell:
'''
##See here for further definition of nested.
## https://www.biostars.org/p/272676/
## https://bedops.readthedocs.io/en/latest/content/reference/set-operations/nested-elements.html

## Note the GFF likely has to be sorted beforehand. 

mkdir filtered

##Basically a reimplmentation of the awk script
##But, had to modify it so it not just "previous feature", but all previous features where prev_end >= cur_start
python -c "
r_handle = open('!{inputGff}')
w_handle = open('nested.gff','w')
prev_features = dict()

SKIP_TYPES=['polypeptide']

nested_by_accept_filter = {'G3DSA:3.40.50.720':'G3DSA:3.90.180.10'} ## Implemented as a heuristic, where only the keys are eligible for nested filtering, since not all nested annotations are so well behaved as these.

for l in r_handle.readlines():
    if l[0] == '#':
        continue
    sl = l.split('\t')
    cur_seqid = sl[0]
    cur_type = sl[2].strip()
    cur_start = int(sl[3])
    cur_end = int(sl[4])
    cur_attr = sl[8].split(';') ## 9th column    
    cur_attr = {x.split('=')[0]:x.split('=')[1] for x in cur_attr}

    if cur_type in SKIP_TYPES:
        print('Skipping',sl)
        continue

    ##Evaluate previous features, drop non relevant ones
    droppable = []
    for f in prev_features:
        if prev_features[f]['end'] < cur_start:
            ##This feature is no longer relevant
            droppable.append(f)
    [prev_features.pop(key) for key in droppable] ##Safely drop all features from the dictionary

    for f in prev_features:
        if (prev_features[f]['start'] <= cur_start) and (prev_features[f]['end'] >= cur_end):
            ##Fully overlapped (nested), left nested, right nested, or fully coincident with 'some feature'
            if cur_attr['Name'] not in nested_by_accept_filter.keys():
                continue
            if nested_by_accept_filter[cur_attr['Name']] == prev_features[f]['attr']['Name']:
                ##Only handle specific enumerated nestings, rather than all possible nestings
                w_handle.write(l) ##Mark for future filtering
                break

    prev_features[l] = {'seqid':cur_seqid,'start':cur_start,'end':cur_end,'attr':cur_attr}
r_handle.close()
w_handle.close()
"

##Remove the nested features.
grep -v -f <(cut -f 1,2,3,4,5 nested.gff) !{inputGff} | grep -vP "^###$" > filtered/!{inputGff}

'''

}

process merge_fasta {
cpus 1
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
input:
 path gff
output:
 path 'ipr_merged.fasta'
shell:
'''
ls -1L | grep ".fasta" > fastas.txt
cat fastas.txt | xargs -n 1 seqkit seq > ipr_merged.fasta
'''
}

process merge_gff {
    cpus 1
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    conda 'genometools-genometools'
    input:
     path collected_Gffs
    output:
     path "ipr_merged.gff"
    shell:
      '''
      ls -1L | grep ".gff" > gffs.txt
      ##cat $(head -n 1 gffs.txt) | grep "##" | grep -v "##FASTA" | grep -v "##sequence-region" > tmp.ipr_merged.gff ##Just the header lines
      ##cat $(head -n 1 gffs.txt) | grep -v "##" >> tmp.ipr_merged.gff ##The content
      cat gffs.txt | xargs -n 1 cat | grep -Pv "^#" >> tmp.ipr_merged.gff
      ##tail -n +2 gffs.txt | xargs -n 1 cat | grep -v "##" >> tmp.ipr_merged.gff
      gt gff3 -tidy -sort -retainids tmp.ipr_merged.gff > ipr_merged.gff
      '''
}

process gt_extractfeat {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
conda 'genometools-genometools seqkit'
input:
 path gff
 path fasta
 val type
output:
 path "${gff}.fa"
shell:
'''
mkdir results
types=("hmmsearch" "protein_match")
for t in ${types[@]}; do
   gt extractfeat -matchdescstart -type ${t} -retainids -coords -seqfile !{fasta} !{gff} > results/!{gff}.${t}.fa
done

##seqkit replace is done, as some programs don't like colons, e.g. raxml-ng
seqkit seq ./results/*.fa | seqkit replace -p ":" -r "=" > !{gff}.fa
'''
}


process ipr_svg {
    cpus 1
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
     path xml
    output:
     path "*.svg"
    shell:
      '''
      ##Note: SVG output from interproscan is deprecated, and no longer mentioned in the documentation
      /oasis/tscc/scratch/tfallon/software/interproscan-5.63-95.0/interproscan.sh -mode convert -d ./ -i !{xml} -f svg
      tar xzf *.tar.gz
      '''
}

process ipr_shorten_gff {
    cpus 1
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
     path gff
    output:
     path "results/${gff}"
    shell:
      '''
      mkdir results
      cat !{gff} | sed -E 's^,"MetaCyc.+^^g' > results/!{gff}
      '''
}


process gt_sketch_svg {
    cpus 1
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
     tuple path(gff),path(style_file)
    output:
     path "*.gt.svg"
    shell:
      '''
      ### UPDATE, unfortunately, I can't get the scale right in gt_sketch, so will try another program.  Namely: https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
      ##Default style file is in /home/tfallon/miniconda3/share/genometools/gtdata/sketch/default.style
      PEP_LENGTH=`cat !{gff} | grep polypeptide | cut -f 5`
      PEP_ID=`cat !{gff} | grep polypeptide | cut -f 1`
      echo "PEP_LENGTH:${PEP_LENGTH}"
      ##Scale the width
      MARGIN=30 ##This width must include the left/right margin in in style_file
      WIDTH=$(python -c "print(int(${PEP_LENGTH} * 2 - $MARGIN * 2))")
      echo "WIDTH:${WIDTH}"
      gt sketch -style !{style_file} -width $WIDTH -format svg ${PEP_ID}.tmp.svg <(cat !{gff} | grep -P "##|SUPERFAMILY")
      ##Below is a workaround since svgutils doesn't rename glyphs
      cat ${PEP_ID}.tmp.svg | sed "s/glyph/glyph_${PEP_ID}/g" > ${PEP_ID}.gt.svg
      '''
}

process extract_non_annotated {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
input:
 path gff
 path fasta
 val minlength
output:
 path('non_annotated-lenfilter.fa')
 path('non_annotated-lenfilter.gff')
conda 'bedops bedtools grep agat seqkit genometools-genometools'
shell:
'''
cat !{gff} | grep "md5=" | gff2bed > polypeptide.bed

bedtools subtract -a polypeptide.bed -b <(cat !{gff} | grep -vP "\tpolypeptide\t") | bedtools sort > non_annotated.bed
cat non_annotated.bed | awk '($3-$2) >=!{minlength}' > non_annotated-lenfilter.bed
agat_convert_bed2gff.pl --bed non_annotated-lenfilter.bed -o non_annotated-lenfilter.tmp.gff

##Just rename gene to polypeptide
cat non_annotated-lenfilter.tmp.gff | sed 's/\tgene\t/\tpolypeptide\t/g' | gt gff3 -tidy -sort -retainids > non_annotated-lenfilter.gff

seqkit subseq --gtf non_annotated-lenfilter.gff !{fasta} > non_annotated-lenfilter.fa

##gt extractfeat -matchdescstart -type 'polypeptide' -retainids -coords -seqfile !{fasta} non_annotated-lenfilter.gff

## The below, is a weird way, to do a coordinate liftover from bed to gff-style coordinates, if outputting a bed instead of a gff
#bedtools getfasta -fi !{fasta} -bed non_annotated-lenfilter.bed > non_annotated-lenfilter.tmp.fa

## A workaround to liftover the bed coordinates from 0-indexed to 1-indexed
## A pox on whoever descided bed should be 0-indexed
#cat non_annotated-lenfilter.tmp.fa | grep ">" | cut -f 2 -d ':' | awk -F'-' -v OFS='-' '{print ($1),$2 "\011" ($1+1),($2+1)}' > mapping.txt
#seqkit replace -p '([0-9]+-[0-9]+)' -r '{kv}' -k mapping.txt non_annotated-lenfilter.tmp.fa > non_annotated-lenfilter.fa
'''
}

process filter_ACP_domain {
input:
 path(fasta)
output:
 path("results/${fasta}")
shell:
'''
##See https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036736/, "Contributing Member Database Entries"
mkdir results
seqkit grep -n -r -p "G3DSA:1.10.1200.10" !{fasta} > results/!{fasta}
'''
}

workflow getTargets {
 take: 
   targets
   fasta 
 main:
 download_pfam_A | prepare_pfam_A
 targetHmms = Channel.fromList(targets)
 hmm_subset(targetHmms,prepare_pfam_A.out.notpressed) 
 hmm_subset.out.collectFile(name: 'targets.hmm') | hmmpress_generic
 detab(fasta) //ran into a FASTA file where the ID was separated with tabs, instead of spaces. It broke seqkit grep -f
 debar(detab.out)
 //decolon(debar.out)
 pfam_run(debar.out.combine(hmmpress_generic.out))
 emit:
      pfam_run.out.matched
}

workflow pfam_AB_run {
take:
   peptideFastas
main:
    pfamA = channel.fromPath("/home/tfallon/databases/pfam/hmmpress/Pfam-A/Pfam-A.hmm.*").toSortedList()
    //pfamB = channel.fromPath("/home/tfallon/databases/pfam/hmmpress/Pfam-B/Pfam-B.hmm.*").toSortedList()
    //pfams = pfamA.mix(pfamB)
    //pfam_run(peptideFastas.combine(pfams))
    pfam_run(peptideFastas.combine(pfamA))
    hmmersearch2gff(pfam_run.out.tbl)
    gff_score_filter(hmmersearch2gff.out,'1E-3')
    //gff_nested_filter(gff_score_filter.out)
emit:
    gff_score_filter.out
}

workflow hmmer_by_pfam {
  fasta = channel.fromPath(params.fasta)
  fastaChunks = fasta.splitFasta(by:700,file:true)
  pfam_AB_run(fastaChunks)
  merge_gff(pfam_AB_run.out.collect())  
  generate_svg_colors(merge_gff.out)
  split_gff_by_seqid(merge_gff.out)
  dummy_publish_path(channel.fromPath('renaming.txt'))
  gff_w_color = split_gff_by_seqid.out.flatten().combine(generate_svg_colors.out).combine(dummy_publish_path.out)
  DNA_features_viewer(gff_w_color)
}

workflow download_data {
 download_interproscan_docker_data()
}

workflow {
    //FASTA should be peptides, ideally already filtered down to those with hits.
    //See the getTargets subworkflow, to automate target searching and selection
    //PF00550 = "PP-binding", ACP overlapping
    //PF14573 = "PP-binding_2", ACP overlapping
    //PF00109 = ketpacyl-synt, KS overlapping
    //PF02801 = Ketoacyl-synt_C, KS overlapping
    //PF16197 = KAsynt_C_assoc, KS overlapping
    //PF00108 = Thiolase_N, KS overlapping
    //PF02803 = Thiolase C, KS overlapping
    //PF01154 = HMG_CoA_synt_N, KS overlappping
    //PF08540 = HMG_CoA_synt_C
    //PF00195 = Chal_sti_synt_N
    //PF02797 = Chal_sti_synt_C

    //download_ipr2go() //Vestigial

    unfiltered = channel.fromPath(params.fasta)
    //getTargets(["PF00550","PF14573","PF00109","PF02801","PF16197","PF00108","PF02803","PF01154","PF08540","PF00195","PF02797"],unfiltered)
    //data = getTargets.out //If running on an unfiltered dataset, this uses HMMER3/PFAM searches to filter out candidates
    data = unfiltered

    fasta_remove_asterisk(data)
    fastaChunks = fasta_remove_asterisk.out.splitFasta(by:700,file:true) 
    //pfam_AB_run(fastaChunks) // Not necesarry for pure interproscan annotation
    interproscan_run(fastaChunks)

    //mgkit_hmmer2gff(pfam_run.out)    

    gff_strip_fasta(interproscan_run.out.gffs)
    gff_append_name_to_seqid_ipr(gff_strip_fasta.out.gffs)
    stripped_gffs = gff_append_name_to_seqid_ipr.out // .mix(pfam_AB_run.out) // Uncomment if wanting to include pfam search results in things.
    collected_gffs = stripped_gffs.collect()
    collected_ipr_fa = gff_strip_fasta.out.fastas.collect()

    collected_jsons = interproscan_run.out.jsons
    tsv_results = interproscan_run.out.tsvs.collectFile()
    sites_results = interproscan_run.out.sites.collectFile()

    //ipr_svg(interproscan_run.out.xmls) // Deprecated, InterProScan no longer supports SVG generation.
    //gt_style = channel.fromPath('ipr.style')
    //gt_sketch_svg(stripped_gffs.combine(gt_style))
    //merge_fasta(collected_ipr_fa)
    
    merge_gff(collected_gffs)

    ipr_shorten_gff(merge_gff.out)

    
    //TODO, evaluate if this gff_nester_filter is actually operating properly. It probably makes at least the sorting assumption 
    gff_nested_filter(ipr_shorten_gff.out)

    generate_svg_colors(gff_nested_filter.out.remaining_gffs)
    split_gff_by_seqid(gff_nested_filter.out.remaining_gffs)

    gff_w_color_renaming = split_gff_by_seqid.out.flatten().combine(generate_svg_colors.out).combine(channel.fromPath(params.renaming))
    
    DNA_features_viewer(gff_w_color_renaming)
    
    gt_extractfeat(gff_nested_filter.out.remaining_gffs,data,"protein_match")

    svg_utils_merge(DNA_features_viewer.out.collect())
    
    svg_resize_page(DNA_features_viewer.out)
    svg_2_pdf(svg_resize_page.out)

    pdf_2_PDF_A_1B(svg_2_pdf.out)
    //extract_non_annotated(merge_gff.out,unfiltered,'90')
}

workflow do_simple_nonparallel_scan {
    //For clarity, and in contrast to workflow 'scan_and_plot', this workflow only does the interproscan annotation, not the plotting
    //It also doesn't have the commented lines to do a HMMER/PFAM target selection.
    //Or otherwise use PFAM.
    //A more elegant way would be to call do_scan from scan_and_plot, but this gets the job done

    //Update: I'm not totally confident that the parallelization is working right. It seems broken when attempted, with GFF ID= or Name= collisions happening. 
    //So, will rename this, and rerun it.

    //download_ipr2go()  //Vestigial

    unfiltered = channel.fromPath(params.fasta)
    data = unfiltered

    fasta_remove_asterisk(data)
    //fastaChunks = fasta_remove_asterisk.out.splitFasta(by:700,file:true) 
    interproscan_run(fasta_remove_asterisk.out)

    gff_strip_fasta(interproscan_run.out.gffs)
    gff_append_name_to_seqid_ipr(gff_strip_fasta.out.gffs)
    stripped_gffs = gff_append_name_to_seqid_ipr.out
    collected_gffs = stripped_gffs.collect()
    collected_ipr_fa = gff_strip_fasta.out.fastas.collect()

    collected_jsons = interproscan_run.out.jsons
    tsv_results = interproscan_run.out.tsvs.collectFile()
    sites_results = interproscan_run.out.sites.collectFile()

    merge_gff(collected_gffs)

    ipr_shorten_gff(merge_gff.out)
    
    gff_nested_filter(ipr_shorten_gff.out)

    emit:
       gff_nested_filter.out.remaining_gffs
}
