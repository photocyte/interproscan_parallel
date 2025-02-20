nextflow.enable.dsl=2

process generate_svg_colors {
publishDir "results",overwrite:true,mode:'link'
conda 'matplotlib'
cache 'deep'
input:
 path gff
output:
 path "feature_to_color.tsv"
shell:
'''
cat !{gff} | grep -oP "Name=[^;]+" | cut -f 2 -d = | sort | uniq > named_features.txt
python -c "
import matplotlib.colors
import os
color_lookup = list(matplotlib.colors.CSS4_COLORS.keys())
i=0
key_val = dict()
r_handle = open('named_features.txt','r')
for l in r_handle.readlines():
    key_val[l.strip()] = color_lookup[i]
    i+=1
    if i >= len(color_lookup):
        i=0
r_handle.close()
w_handle = open('feature_to_color.tsv','w')
for k in key_val:
    w_handle.write('\\t'.join([k,key_val[k]])+os.linesep)
w_handle.close()
"
'''
}

process svg_2_pdf {
//Running off a singularity container
publishDir "results/${task.process}",overwrite:true,mode:'link'
input:
 path svg
output:
 path "${svg}.pdf"
//From https://superuser.com/questions/381125/how-do-i-convert-an-svg-to-a-pdf-on-linux
script:
"""
##rsvg-convert -f pdf -o ${svg}.pdf ${svg}
inkscape ${svg} --export-pdf=${svg}.pdf
"""
}

process pdf_2_PDF_A_1B {
conda 'ghostscript'
publishDir "results/${task.process}",overwrite:true,mode:'link'
input:
 path pdf
output:
 path "output/${pdf}"
script:
"""
wget https://raw.githubusercontent.com/matteosecli/pdf2archive/master/pdf2archive
chmod +x pdf2archive
mkdir -p output
./pdf2archive ${pdf} output/${pdf}
"""
}

process svg_utils_merge {
 publishDir "results/${task.process}",overwrite:true,mode:'link'
 conda 'svgutils' 
 input:
  path collected_SVGs
 output:
  path "*_compose.svg"
shell:
'''
#! /usr/bin/env python
import svgutils
import svgutils.compose
##import svgutils.transform
import glob

files = sorted(glob.glob("*.svg"))

merge_name = "_".join(files)

##Note the *[] syntax,
##https://stackoverflow.com/questions/3941517/converting-list-to-args-when-calling-function
svgutils.compose.Figure("30cm", str(len(files)*1)+"cm",
        *[ svgutils.compose.Panel(svgutils.compose.SVG(f)) for f in files ]
        ).tile(1, len(files)).save(merge_name+"_compose.svg")
### ^^ Above "compose" approach seems to scramble the text?
## UPDATE: above problem was due to svgutils not renaming the glyphs
## Seems to work okay with the matplotlib derived svgs

#import svgutils.transform
#import glob
#files = glob.glob("*.svg")
##files = files[0:2]

#fig = svgutils.transform.SVGFigure("100cm", str(len(files)*3)+"cm")

#panels = [ svgutils.transform.fromfile(f) for f in files ]

#plots = []
#for i in range(0,len(panels)):
#    plots.append(panels[i].getroot())
#    plots[-1].moveto(0,i*100,scale=1.0) ##moveto works in place.
     
#fig.append(plots)
#fig.save("fig_final.svg")

'''
}

process svg_resize_page {
//conda "inkscape"
//Runs off a container instead.
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
input:
 path svg
output:
 path "mod/$svg"
shell:
'''
mkdir mod
inkscape --actions="select-all; transform-scale:0.01; fit-canvas-to-selection; export-filename:mod/!{svg}; export-do;" !{svg}

##inkscape --export-type="svg" --export-area-drawing -o mod/!{svg} !{svg}
'''
}

process split_gff_by_seqid {
    cpus 1
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    conda 'genometools-genometools'
    cache 'deep'
    input:
     path merged_gff
    output:
     path "split/*.gff" //Can't do gff3, as DNA_Features_Viewer doesn't recognize
    shell:
      '''
      mkdir split
      cat !{merged_gff} | grep -v "#" | cut -f 1 | sort | uniq > seqids.txt
      while read s; do
      cat !{merged_gff} | grep $s | gt gff3 -tidy -sort -retainids > split/$s.gff
      done <seqids.txt
      '''
}

process centroid_gff_features {
publishDir "results/${task.process}", mode: 'link',overwrite:'true'
input:
 path input_gff
output:
 path "mod/${input_gff}"
shell:
'''
mkdir mod
cat !{input_gff} | awk -v "FS=\t" -v "OFS=\t" '{print $1,$2,$3,int(($4+$5)/2-30),int(($4+$5)/2)+30,$6,$7,$8,$9}' > mod/!{input_gff}
'''
}

process equally_distribute_gff_features {
publishDir "results/${task.process}", mode: 'link',overwrite:'true'
conda 'numpy'
input:
 path input_gff
output:
 path "mod/${input_gff}"
shell:
'''
#! /usr/bin/env python
import numpy
import os
starts = []
stops = []

os.mkdir('mod')

r_handle = open("!{input_gff}","r")
w_handle = open("mod/!{input_gff}","w")

spacing = 65

for l in r_handle.readlines():
    split_l = l.split("\\t")
    if l[0] == "#":
        continue
    starts.append(int(split_l[3]))
    stops.append(int(split_l[4]))
r_handle.close()
the_linspace = numpy.linspace(starts[0],starts[-1],num=len(starts))
norm_linspace= ( the_linspace- min(the_linspace) )/( max(the_linspace)-min(the_linspace) )*len(starts)*(spacing+9) ## https://www.geeksforgeeks.org/how-to-normalize-an-numpy-array-so-the-values-range-exactly-between-0-and-1/
norm_starts=norm_linspace.astype(int)
norm_ends = (norm_linspace+spacing).astype(int)
norm_pairs = list(zip(norm_starts,norm_ends))

r_handle = open("!{input_gff}","r")
i=0
for l in r_handle.readlines():
    split_l = l.split("\\t")
    if l[0] == "#":
        continue
    w_handle.write("\t".join([split_l[0],split_l[1],split_l[2],str(norm_pairs[i][0]),str(norm_pairs[i][1]),split_l[5],split_l[6],split_l[7],split_l[8]])+os.linesep) 
    i+=1
w_handle.close()
'''
}

process split_and_restart_gff_by_max_features {
publishDir "results/${task.process}", mode: 'link',overwrite:'true'
input:
 path input_gff
output:
 path "mod/*"
shell:  
'''
#! /usr/bin/env python
import os
os.mkdir('mod')

max_feature_num = 30
subtract_offset = 0
r_handle = open("!{input_gff}","r")
i=1
fnum=1
w_handle = open('mod/{index}_!{input_gff}'.format(index=fnum),'w')
for l in r_handle.readlines():
    split_l = l.split("\\t")
    if l[0] == "#":
        continue
    w_handle.write("\t".join([split_l[0],split_l[1],split_l[2],str(int(split_l[3])-subtract_offset),str(int(split_l[4])-subtract_offset),split_l[5],split_l[6],split_l[7],split_l[8]])+os.linesep)
    if i == max_feature_num:
        w_handle.close()
        fnum += 1
        subtract_offset = int(split_l[4])+9
        w_handle = open('mod/{index}_!{input_gff}'.format(index=fnum),'w')
        i=0
    i+=1 
'''
}

process DNA_features_viewer {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    conda 'matplotlib dna_features_viewer bcbio-gff'
    cache 'deep'
    input:
     tuple path(gff),path(colors),path(renaming)
    output:
     path("*.dfv.svg")
    tag "${gff}"
shell:
'''
#! /usr/bin/env python
##installed with:
##pip install dna_features_viewer
##pip install bcbio-gff
## Note the additional backslash for Nextflow escaping

doRename = True ## Used for debugging of the renaming process. False=Print the method identifier for each technique.

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
plt.rcParams['svg.fonttype'] = 'none'

pep_name = None

global accept_all
accept_all = False
global accepted_sources
accepted_sources = ["GENE3D","SUPERFAMILY","PFAM","PROSITEPATTERNS","PROSITEPROFILES","PIRSR","CDD","SMART"] ## Only used for naming of the domains. They should all graphically be there!
global accepted_types
accepted_types = ["protein_match","hmmsearch"] ## Only used for coloration of the domains. They should all be there!

##This just sets the pep_name from the seqid in the gff file
r_handle = open("!{gff}","r")
for l in r_handle.readlines():
    split_l = l.split("\\t")
    if l[0] == "#":
        continue
    if len(split_l) > 2 and (split_l[2] in accepted_types or accept_all):
        pep_name = split_l[0]
    if split_l[2] == "polypeptide":
        seq_len = int(split_l[4]) - int(split_l[3])
r_handle.close()

##If the pep_name wasn't able to be set from the gff lines, just use the filename
if pep_name == None or pep_name == "":
    pep_name = "!{gff}"
    pep_name = pep_name[:-4] ##Strip off the .gff suffix

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
        if (feature.type in accepted_types and (feature.qualifiers["source"][0].upper() in accepted_sources)) or accept_all:
             return color_kv[feature.qualifiers["Name"][0]].strip()
        else:
            return "gold" ##Just the color

    def compute_feature_label(self, feature):
        if feature.type == 'polypeptide':
            return feature.qualifiers["ID"][0]
        elif (feature.type in accepted_types and (feature.qualifiers["source"][0].upper() in accepted_sources)) or accept_all:
            theName = feature.qualifiers["Name"][0]
            if doRename and theName in rename_kv.keys():
                return rename_kv[theName]
            else:
                return feature.qualifiers["Name"][0] ##Can also do "Name"
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

graphic_record = IPSCustomTranslator().translate_record("!{gff}")
print(f'graphic_record.sequence:{graphic_record.sequence}')
if graphic_record.sequence == None:
    plot_len = seq_len/100.0
else:
    plot_len = len(graphic_record.sequence)/100.0
print(dir(graphic_record))
ax, _ = graphic_record.plot(figure_width=plot_len, strand_in_label_threshold=7)
##ax.figure.savefig(pep_name+'.dfv.svg', bbox_inches='tight') # SAVE AS SVG
ax.figure.savefig('!{gff}.dfv.svg', bbox_inches='tight', transparent=True) # SAVE AS SVG
'''
}

workflow beads_on_string {
    gff_ch = Channel.fromPath(params.gff)
    renaming_ch = Channel.fromPath(params.renaming)
    colors_ch = Channel.fromPath(params.colors)
    centroid_gff_features(gff_ch)
    equally_distribute_gff_features(centroid_gff_features.out)

    //generate_svg_colors(equally_distribute_gff_features.out)
    split_gff_by_seqid(equally_distribute_gff_features.out)

    split_and_restart_gff_by_max_features(split_gff_by_seqid.out)
    
    gff_w_color = split_and_restart_gff_by_max_features.out.flatten().combine(colors_ch).combine(renaming_ch)
    DNA_features_viewer(gff_w_color)

    svg_utils_merge(DNA_features_viewer.out.collect())
    svg_resize_page(svg_utils_merge.out)
    svg_2_pdf(DNA_features_viewer.out.mix(svg_resize_page.out))
    pdf_2_PDF_A_1B(svg_2_pdf.out)

}

workflow {
    gff_ch = Channel.fromPath(params.gff)
    centroid_gff_features(gff_ch)
    generate_svg_colors(gff_ch)
    split_gff_by_seqid(gff_ch)
    gff_w_color = split_gff_by_seqid.out.flatten().combine(generate_svg_colors.out).combine(channel.fromPath('renaming.txt'))
    DNA_features_viewer(gff_w_color)

    svg_utils_merge(DNA_features_viewer.out.collect())
    svg_2_pdf(DNA_features_viewer.out)
    pdf_2_PDF_A_1B(svg_2_pdf.out)
}
