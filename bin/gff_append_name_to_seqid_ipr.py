#!/usr/bin/env python
import re
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("gff")
args = parser.parse_args()

if not os.path.exists('results'):
    os.mkdir('results')

##Below targets no longer used. It is assumed everything will be renamed.
##Check git if wanting to revert the behavior
##targets = ["Gene3D","SUPERFAMILY"] ## Where the rename targets are set. Others are left untouched.

r_handle = open(args.gff)
w_handle = open(f'results/{args.gff}','w')
for l in r_handle.readlines():
    if l[0] == "#":
        w_handle.write(l)
        continue
    NAME_MATCH = re.search("Name=([^;$]+)",l)
    if NAME_MATCH != None:
        NAME=NAME_MATCH.group(1)
    else:
        w_handle.write(l)
        continue
    SEQID = l.split("\\t")[0]
    new_l = re.sub('match\\$[0-9]+_',SEQID+"__"+NAME+"_",l) ## Unclear what the numbers after match$ mean, so will delete.
    w_handle.write(new_l)
r_handle.close()
w_handle.close()
