#!/usr/bin/env python
"""
usage:

Rename prodigal FASTA sequences by adding sample name.
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-p", "--genesPath", required=True)
ap.add_argument("-n", "--sampleName", required=True)
args = vars(ap.parse_args())

with open(args['genesPath'], 'r') as f:
    lines = []
    for line in f:
        if line.startswith(">"):
            split_line = line.split(maxsplit=1)
            split_line[0] = split_line[0] + '_' + args['sampleName']
            lines.append(' '.join(split_line))
        else:
            lines.append(line)

file_name = args['genesPath'].split(sep='/')[-1].split(sep='.')
file_name[0] = file_name[0] + '_' + args['sampleName']
file_name = '.'.join(file_name)

with open(file_name, 'w') as f:
    for line in lines:
        f.write(line)
