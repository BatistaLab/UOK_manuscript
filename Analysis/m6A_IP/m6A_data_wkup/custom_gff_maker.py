#!/usr/bin/env python3
# This script adds quotes to the attributes column to generate a custom gff file

from pathlib import Path

p = Path('input.txt')

def convert_line(line):
    fields_in = line.strip().split('\t')
    id_in = fields_in[8]

    id_out = f'gene_id "{id_in}"'
    fields_out = [*fields_in[:-1], id_out]

    return '\t'.join(fields_out) + '\n'

f_in = open('input.txt')
f_out = open('out.gff', 'w')

for line_in in f_in:
    line_out = convert_line(line_in)
    f_out.write(line_out)
