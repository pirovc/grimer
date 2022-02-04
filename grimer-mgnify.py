#!/usr/bin/env python3

import scripts.mgnify_download
import grimer.grimer
import argparse
import os
import glob

parser = argparse.ArgumentParser(description='grimer-mgnify')
parser.add_argument('-i', '--mgnify-study-accession', required=True, type=str, help="MGnify study accession (e.g. MGYS00002462)")
parser.add_argument('-g', '--grimer-params', type=str, help="Extra params for grimer")
parser.add_argument('-o', '--output-prefix', type=str, help="Output prefix for files and report")
args = parser.parse_args()

if args.output_prefix:
    prefix = args.output_prefix
else:
    prefix = args.mgnify_study_accession

# download files
print("Downloading files for study accession " + args.mgnify_study_accession)
scripts.mgnify_download.main(['-i', args.mgnify_study_accession, '-o', prefix, '-v'])

files = filter(os.path.isfile, glob.glob(prefix + '*taxonomy_abundances*'))
# Sort files by size ASC
files = sorted(files, key=lambda x: os.stat(x).st_size)
md = glob.glob(prefix + '*_metadata.tsv*')

grimer.grimer.main(["-i", files[-1],
                    "-m", md[-1],
                    "-c", 'config/default.yaml',
                    "-f", ";",
                    "--obs-replace", "^.+__", "", "_", " ",
                    "-r", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species",
                    "-t", "ncbi",
                    "-o", prefix + ".html",
                    "--title", "MGnify study accession " + args.mgnify_study_accession,
                    *args.grimer_params.split(" ")
                    ])
