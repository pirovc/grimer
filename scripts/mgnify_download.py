#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import pickle
import gzip
from urllib.parse import urlencode
from jsonapi_client import Session, Filter
from glob import glob

"""
Script to download taxonomy abundance files and metadata from MGnify
- It always downloads latest available results based on the pipeline version

Usage: ./mgnify_download.py study_accession [output_folder]
Example single: ./mgnify_download.py MGYS00000554
Example dump: seq -f "MGYS%08g" 1 5724 | xargs -P 8 -I {} ./mgnify_download.py -i {} -v -g -o mgnify_dump_20210408/ > mgnify_dump_20210408.log 2>&1 &
"""


def main(argv=sys.argv[1:]):

    API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/'

    parser = argparse.ArgumentParser(description='grimer-download-mgnify')
    parser.add_argument('-i', '--study-accession', required=True, type=str, help="MGnify study accession (e.g. MGYS00002462)")
    parser.add_argument('-g', '--gzip', default=False, action='store_true', help="Gzip downloaded files")
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help="Verbose output")
    parser.add_argument('-o', '--output-prefix', type=str, help="Output prefix for downloaded files. Default: --study-accession")
    args = parser.parse_args(argv)

    study_accession = args.study_accession
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        prefix = study_accession
    gz = args.gzip

    md_file = prefix + "_metadata.tsv"
    out_file = prefix + ".pkl"
    if gz:
        out_file = out_file + ".gz"
        md_file = md_file + ".gz"

    # Check if files exists and skip
    tax_files = glob(prefix + "*_taxonomy_abundances_*")
    if tax_files and os.path.isfile(out_file) and os.path.isfile(md_file):
        print(study_accession, "Warning: files already exist ")
        return

    with Session(API_BASE) as s:
        # Get main study resource
        try:
            study = s.get('studies', study_accession).resource
            if args.verbose:
                print(study.accession, "SAMPLES:" + str(study.samples_count), sep="\t", end="\t")
        except:
            print(study_accession, "Error: Study accession not found")
            sys.exit(1)

        # Save study info as a dict in a pkl file
        f = gzip.open(out_file, 'wb') if gz else open(out_file, "wb")
        pickle.dump(study.json, file=f)
        f.close()

        # Get all taxonomic tables for the highest version of the pipeline
        highest_version = 0
        table_version = {}
        for download in study.downloads:
            label = download.description.label
            #["Taxonomic assignments",
            #"Taxonomic assignments SSU",
            #"Taxonomic assignments LSU"
            #"Taxonomic assignments UNITE",
            #"Taxonomic assignments ITSoneDB"]
            if "Taxonomic assignments" in label:
                version = float(download.pipeline.id)
                if version not in table_version:
                    table_version[version] = []
                table_version[version].append(download.url)
                if version > highest_version:
                    highest_version = version

        if not table_version:
            print("Error: No taxonomic assignments for this study to download")
            sys.exit(1)
        else:
            table_urls = table_version[highest_version]

        # Get all available samples in one go and collect metadata
        params = {
            'study_accession': study_accession,
            'page_size': study.samples_count,
        }
        fltr = Filter(urlencode(params))

        metadata = {}
        for sample in s.iterate('samples', fltr):
            # TODO: how to access runs faster, sample.runs is too slow
            #nruns += len(sample.runs)
            metadata[sample.accession] = {}
            for md in sample.sample_metadata:
                metadata[sample.accession][md["key"]] = md["value"]
            # Add sample description, name and name as metadata
            metadata[sample.accession]['sample-desc'] = sample.sample_desc
            metadata[sample.accession]['sample-name'] = sample.sample_name

        # Get link sample accession, run accession
        # TODO treat multiple runs per sample
        run_sample_accesion = {}
        try:
            for run in s.iterate('runs', fltr):
                run_sample_accesion[run.sample.id] = run.id
        except:
            print("Error: Could not retrieve run accession", sep="\t", end="\t")

    # Write metadata
    md_df = pd.DataFrame.from_dict(metadata).T
    if run_sample_accesion:
        mapped_accessions = md_df.index.isin(run_sample_accesion.keys())
        if args.verbose:
            print("MAPPED:" + str(sum(mapped_accessions)), sep="\t", end="\t")
        md_df.index = md_df.index.map(lambda x: run_sample_accesion[x] if x in run_sample_accesion else x)
    else:
        if args.verbose:
            print("Warning: No mapping between accessions of samples and metadata", sep="\t", end="\t")

    if args.verbose:
        print("METADATA:" + str(md_df.shape[1]), sep="\t", end="\t")
    md_df.to_csv(md_file, compression="gzip" if gz else None, sep="\t")

    # Read and write tables
    for table_url in table_urls:
        try:
            t = pd.read_table(table_url)
            if args.verbose:
                print("OK:" + table_url, end=";")
            # Print original
            filename = prefix + "_" + os.path.basename(table_url)
            t.to_csv(filename if not gz else filename + ".gz", compression="gzip" if gz else None, sep="\t", index=False)
        except:
            if args.verbose:
                print("INVALID:" + table_url, end=";")

    if args.verbose:
        print()


if __name__ == "__main__":
    main()
