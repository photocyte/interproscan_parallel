// Also, see here https://lipm-gitlab.toulouse.inra.fr/LIPM-BIOINFO/nextflow-interproscan5/blob/master/interproscan.nf
nextflow.enable.dsl=2

include { download_pfam_A } from '/home/tfallon/databases/pfam/pfam.nf'
include { prepare_pfam_A } from '/home/tfallon/databases/pfam/pfam.nf'
include { pfam_run } from '/home/tfallon/databases/pfam/pfam.nf'
include { hmmpress_generic } from '/home/tfallon/databases/pfam/pfam.nf'
include { generate_svg_colors ; svg_2_pdf ; svg_utils_merge ; split_gff_by_seqid ; DNA_features_viewer } from '/home/tfallon/source/nextflow/interproscan_parallel/dna_features_viewer.nf'

process interproscan_run {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    //storeDir "results/interproscan_run"
    cpus 1
    //scratch true
    scratch 'ram-disk'
    //stageInMode 'copy'
    errorStrategy 'finish'
    maxRetries 3
    executor "local"
    queue "home"
    time '24h'
    maxForks 10
    cache "deep"
    input:
     path inputFasta
    output:
     path "*.gff3", emit: gffs
     path "*.json", emit: jsons
     path "*.tsv", emit: tsvs 
     path "*.tsv.sites", emit: sites, optional:true
     path "*.xml", emit: xmls
    script:
      """
      /home/tfallon/software/interproscan-5.51-85.0/interproscan.sh --disable-precalc --iprlookup --pathways --goterms --enable-tsv-residue-annot --cpu ${task.cpus} -i ${inputFasta} --tempdir ./
      """
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
##Note the additional backslash for Nextflow escaping
mkdir results
seqkit replace -s -p "\\*$" -r "" !{fasta} > results/!{fasta}
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
TARGETS="##|Gene3D|polypeptide|SUPERFAMILY"
head -n +`grep -n "##FASTA" !{gff} | cut -d : -f 1` !{gff} | grep -v "##FASTA" | grep -P ${TARGETS} > results/ipr.!{gff}
tail -n +`grep -n "##FASTA" !{gff} | cut -d : -f 1` !{gff} | grep -v "##FASTA" > results/ipr.!{gff}.fasta
'''
}

process gff_rename_ipr {
cpus 1
input:
 path gff
output:
 path "results/${gff}"
shell:
'''
#!/usr/bin/env python
import re
import os
os.mkdir('results')
targets = ["Gene3D","SUPERFAMILY"]
r_handle = open('!{gff}')
w_handle = open('results/!{gff}','w')
for l in r_handle.readlines():
    if l[0] == "#":
        w_handle.write(l)
        continue
    NOMATCH = True
    for t in targets:
        if t in l:
            NOMATCH = False
            NAME = re.search("Name=([\\w:.]+)",l).group(1)
            SEQID = l.split("\\t")[0]
            new_l = re.sub('match\\$',SEQID+"_"+NAME+"_",l)
            w_handle.write(new_l)
    if NOMATCH == True:
        w_handle.write(l)
r_handle.close()
w_handle.close()
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
output:
 path "filtered/${inputGff}"
shell:
'''
mkdir filtered
cat !{inputGff} | bioawk -c gff '$score < 1E-3{ print $0 }' | gt gff3 -tidy -sort -retainids > filtered/!{inputGff}
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
input:
 path inputGff
output:
 path "filtered/${inputGff}"
shell:
'''
##See here for this trick
## https://www.biostars.org/p/272676/
mkdir filtered
cat !{inputGff} | gff2bed > !{inputGff}.bed

##Basically a reimplmentation of the awk script
python -c "
r_handle = open('!{inputGff}.bed')
w_handle = open('nested.bed','w')
prev_start = 0
prev_end = 0
for l in r_handle.readlines():
    if l[0] == '#':
        continue
    sl = l.split('\t')
    cur_start = int(sl[1])
    cur_end = int(sl[2])
    if (prev_start <= cur_start) and (prev_end >= cur_end):
        ##Fully overlapped
        w_handle.write(l)
    prev_start = cur_start
    prev_end = cur_end
r_handle.close()
w_handle.close()
"
cat nested.bed | bedops --not-element-of 100% !{inputGff}.bed - > !{inputGff}.fil.bed
grep -f <(cut -f 4 !{inputGff}.fil.bed) !{inputGff} | gt gff3 -tidy -sort -retainids > filtered/!{inputGff}
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
seqkit seq ./results/*.fa > !{gff}.fa
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
      /oasis/tscc/scratch/tfallon/software/interproscan-5.51-85.0/interproscan.sh -mode convert -d ./ -i !{xml} -f svg
      tar xzf *.tar.gz
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

process DNA_features_viewer {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
     tuple path(gff),path(colors),path(renaming)
    output:
     path "*.dfv.svg"
    tag "${gff}"
    shell:
'''
#! /usr/bin/env python
##installed with:
##pip install dna_features_viewer
##pip install bcbio-gff
## Note the additional backslash for Nextflow escaping

import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
pep_name = None
global accepted_sources
accepted_sources = ["Gene3D","SUPERFAMILY","PFAM"]
global accepted_types
accepted_types = ["protein_match","hmmsearch"]

##This just sets the pep_name from the seqid in the gff file
r_handle = open("!{gff}","r")
for l in r_handle.readlines():
    split_l = l.split("\\t")
    if l[0] == "#":
        continue
    if len(split_l) > 2 and split_l[2] in accepted_types:
        pep_name = split_l[0]
r_handle.close()

##If the pep_name wasn't able to be set from the gff lines, just use the filename
if pep_name == None:
    pep_name = "!{gff}"
    pep_name = pep_name[-3:] ##Strip off the suffix

r_handle = open("!{colors}","r")
color_kv = dict()
for l in r_handle.readlines():
    split_l = l.split("\\t")
    color_kv[split_l[0]] = split_l[1]
r_handle.close()

r_handle = open("!{renaming}","r")
rename_kv = dict()
for l in r_handle.readlines():
    split_l = l.split("\\t")
    rename_kv[split_l[0]] = split_l[1].strip()
r_handle.close()

from dna_features_viewer import BiopythonTranslator

class IPSCustomTranslator(BiopythonTranslator):

    def compute_feature_color(self,feature):
        print(feature)
        print(dir(feature))
        if feature.type in accepted_types and (feature.qualifiers["source"][0] in accepted_sources):
             return color_kv[feature.qualifiers["Name"][0]].strip()
        else:
            return "gold"

    def compute_feature_label(self, feature):
        if feature.type == 'polypeptide':
            return feature.qualifiers["ID"][0]
        elif feature.type in accepted_types and (feature.qualifiers["source"][0] in accepted_sources):
            theName = feature.qualifiers["Name"][0]
            if theName in rename_kv.keys():
                return rename_kv[theName]
            else:
                return feature.qualifiers["Name"][0]
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

graphic_record = IPSCustomTranslator().translate_record("!{gff}")
print(dir(graphic_record))
ax, _ = graphic_record.plot(figure_width=len(graphic_record.sequence)/100.0, strand_in_label_threshold=7)
ax.figure.savefig(pep_name+'.dfv.svg', bbox_inches='tight') # SAVE AS SVG
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
 pfam_run(fasta.combine(hmmpress_generic.out))
 emit:
      pfam_run.out.matched
}

workflow pfam_AB_run {
take:
   peptideFastas
main:
    pfamA = channel.fromPath("/home/tfallon/databases/pfam/hmmpress/Pfam-A/Pfam-A.hmm.*").toSortedList()
    pfamB = channel.fromPath("/home/tfallon/databases/pfam/hmmpress/Pfam-B/Pfam-B.hmm.*").toSortedList()
    pfams = pfamA.mix(pfamB)
    pfam_run(peptideFastas.combine(pfams))
    hmmersearch2gff(pfam_run.out.tbl)
    gff_score_filter(hmmersearch2gff.out)
    gff_nested_filter(gff_score_filter.out)
emit:
    gff_nested_filter.out
}

workflow hmmer_by_pfam {
  fasta = channel.fromPath(params.fasta)
  fastaChunks = fasta.splitFasta(by:10,file:true)
  pfam_AB_run(fastaChunks)
  merge_gff(pfam_AB_run.out.collect())  
  generate_svg_colors(merge_gff.out)
  split_gff_by_seqid(merge_gff.out)
  gff_w_color = split_gff_by_seqid.out.flatten().combine(generate_svg_colors.out).combine(channel.fromPath('renaming.txt'))
  DNA_features_viewer(gff_w_color)
}

workflow {
    //FASTA should be peptides, ideally already filtered down to those with hits.
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
    unfiltered = channel.fromPath(params.fasta)
    getTargets(["PF00550","PF14573","PF00109","PF02801","PF16197","PF00108","PF02803","PF01154","PF08540","PF00195","PF02797"],unfiltered)
    data = getTargets.out

    fasta_remove_asterisk(data)
    fastaChunks = fasta_remove_asterisk.out.splitFasta(by:5,file:true)
    pfam_AB_run(fastaChunks)
    interproscan_run(fastaChunks)

    //mgkit_hmmer2gff(pfam_run.out)    

    gff_strip_fasta(interproscan_run.out.gffs)
    gff_rename_ipr(gff_strip_fasta.out.gffs)
    stripped_gffs = gff_rename_ipr.out.mix(pfam_AB_run.out)
    collected_gffs = stripped_gffs.collect()
    collected_ipr_fa = gff_strip_fasta.out.fastas.collect()

    collected_jsons = interproscan_run.out.jsons
    tsv_results = interproscan_run.out.tsvs.collectFile()
    sites_results = interproscan_run.out.sites.collectFile()

    //ipr_svg(interproscan_run.out.xmls)
    gt_style = channel.fromPath('ipr.style')
    //gt_sketch_svg(stripped_gffs.combine(gt_style))
    //merge_fasta(collected_ipr_fa)
    merge_gff(collected_gffs)
 
    generate_svg_colors(merge_gff.out)
    split_gff_by_seqid(merge_gff.out)
    gff_w_color = split_gff_by_seqid.out.flatten().combine(generate_svg_colors.out).combine(channel.fromPath('renaming.txt'))
    DNA_features_viewer(gff_w_color)
    

    svg_utils_merge(DNA_features_viewer.out.collect())
    svg_2_pdf(DNA_features_viewer.out)
    //gt_extractfeat(merge_gff.out,data,"protein_match")
}
