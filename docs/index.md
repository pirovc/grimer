# GRIMER

<img src="https://raw.githubusercontent.com/pirovc/grimer/main/grimer/img/logo.png">

## About

GRIMER is a tool that performs automated analyses and generates a portable and interactive dashboard integrating annotation, taxonomy and metadata. It unifies several sources of evidence to help detect contamination. GRIMER is independent of quantification methods and directly analyses contingency tables to create an interactive and offline report. Reports can be created in seconds and are accessible for non-specialists, providing an intuitive set of charts to explore data distribution among observations and samples and its connections with external sources.

- More information about the method: [pre-print](https://doi.org/10.1101/2021.06.22.449360){ target="_blank" }
- Source-code: [GitHub repository](https://github.com/pirovc/grimer){ target="_blank" }

## Installation

Via conda

```bash
conda install -c bioconda -c conda-forge grimer
```

or locally installing only dependencies via conda:

```bash
git clone https://github.com/pirovc/grimer.git
cd grimer
conda env create -f env.yaml # or mamba env create -f env.yaml
conda activate grimer # or source activate grimer
python setup.py install --record files.txt # Uninstall: xargs rm -rf < files.txt
grimer -h
```

## Basic Usage

- In-depth examples of input files: [Importing files](importing)
- Complete examples of usage with real files: [Examples](examples)


Tab-separated input table

```bash
grimer -i input_table.tsv
```

BIOM file
```bash
grimer -i myfile.biom
```

Tab-separated input table with taxonomic annotated observations (e.g. sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria...)
```bash
grimer -i input_table.tsv -f ";"
```

Tab-separated input table with metadata
```bash
grimer -i input_table.tsv -m metadata.tsv
```

With taxonomy integration (ncbi)
```bash
grimer -i input_table.tsv -m metadata.tsv -t ncbi #optional -b taxdump.tar.gz
```

With configuration file to setup external tools, references and annotations
```bash
grimer -i input_table.tsv -m metadata.tsv -t ncbi -c config/default.yaml -d -g
```

## Parameters


           ▄████  ██▀███   ██▓ ███▄ ▄███▓▓█████  ██▀███  
          ██▒ ▀█▒▓██ ▒ ██▒▓██▒▓██▒▀█▀ ██▒▓█   ▀ ▓██ ▒ ██▒
         ▒██░▄▄▄░▓██ ░▄█ ▒▒██▒▓██    ▓██░▒███   ▓██ ░▄█ ▒
         ░▓█  ██▓▒██▀▀█▄  ░██░▒██    ▒██ ▒▓█  ▄ ▒██▀▀█▄  
         ░▒▓███▀▒░██▓ ▒██▒░██░▒██▒   ░██▒░▒████▒░██▓ ▒██▒
          ░▒   ▒ ░ ▒▓ ░▒▓░░▓  ░ ▒░   ░  ░░░ ▒░ ░░ ▒▓ ░▒▓░
           ░   ░   ░▒ ░ ▒░ ▒ ░░  ░      ░ ░ ░  ░  ░▒ ░ ▒░
         ░ ░   ░   ░░   ░  ▒ ░░      ░      ░     ░░   ░ 
               ░    ░      ░         ░      ░  ░   ░     
         version 1.1.0


    usage: grimer [-h] -i INPUT_FILE [-m METADATA_FILE] [-c CONFIG]
                  [-t {ncbi,gtdb,silva,greengenes,ott}] [-b [TAXONOMY_FILES ...]] [-r [RANKS ...]]
                  [-l TITLE] [-p [{overview,samples,heatmap,correlation} ...]] [-o OUTPUT_HTML]
                  [--full-offline] [-g] [-d] [-f LEVEL_SEPARATOR] [-y VALUES] [-w] [-s]
                  [-u [UNASSIGNED_HEADER ...]] [-z REPLACE_ZEROS] [--obs-replace [OBS_REPLACE ...]]
                  [--sample-replace [SAMPLE_REPLACE ...]] [--min-frequency MIN_FREQUENCY]
                  [--max-frequency MAX_FREQUENCY] [--min-count MIN_COUNT] [--max-count MAX_COUNT]
                  [-j TOP_OBS_BARS] [-a {none,norm,log,clr}] [-e METADATA_COLS] [--optimal-ordering]
                  [--show-zeros]
                  [--linkage-methods [{single,complete,average,centroid,median,ward,weighted} ...]]
                  [--linkage-metrics [{braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,kulczynski1,mahalanobis,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,yule} ...]]
                  [--skip-dendrogram] [-x TOP_OBS_CORR] [-v]

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

    required arguments:
      -i INPUT_FILE, --input-file INPUT_FILE
                            Tab-separatad file with table with counts (Observation table, Count table,
                            Contingency Tables, ...) or .biom file. By default rows contain observations
                            and columns contain samples (use --transpose if your file is reversed). The
                            first column and first row are used as headers. (default: None)

    main arguments:
      -m METADATA_FILE, --metadata-file METADATA_FILE
                            Tab-separated file with metadata. Rows should contain samples and columns
                            the metadata fields. QIIME2 metadata format is accepted, with an extra row
                            to define categorical and numerical fields. If --input-file is a .biom file,
                            metadata will be extracted from it if available. (default: None)
      -c CONFIG, --config CONFIG
                            Configuration file with definitions of references, controls and external
                            tools. (default: None)
      -t {ncbi,gtdb,silva,greengenes,ott}, --taxonomy {ncbi,gtdb,silva,greengenes,ott}
                            Enable taxonomic analysis, convert entries and annotate samples. Files will
                            be automatically downloaded and parsed. Optionally, stored files can be
                            provided with --taxonomy-files. (default: None)
      -b [TAXONOMY_FILES ...], --taxonomy-files [TAXONOMY_FILES ...]
                            Specific taxonomy files to use with --taxonomy. (default: [])
      -r [RANKS ...], --ranks [RANKS ...]
                            Taxonomic ranks to generate visualizations. Use 'default' to use entries
                            from the table directly. (default: ['default'])

    output arguments:
      -l TITLE, --title TITLE
                            Title to display on the top of the report. (default: )
      -p [{overview,samples,heatmap,correlation} ...], --output-plots [{overview,samples,heatmap,correlation} ...]
                            Plots to generate. (default: ['overview', 'samples', 'heatmap',
                            'correlation'])
      -o OUTPUT_HTML, --output-html OUTPUT_HTML
                            Filename of the HTML report output. (default: output.html)
      --full-offline        Embed Bokeh javascript library in the output file. Output will be around
                            1.5MB bigger but it will work without internet connection. ~your report will
                            live forever~ (default: False)

    general data options:
      -g, --mgnify          Plot MGnify, requires --config file with parsed MGnify database. (default:
                            False)
      -d, --decontam        Run DECONTAM and generate plots. requires --config file with DECONTAM
                            configuration. (default: False)
      -f LEVEL_SEPARATOR, --level-separator LEVEL_SEPARATOR
                            If provided, consider --input-table to be a hierarchical multi-level table
                            where the observations headers are separated by the indicated separator char
                            (usually ';' or '|') (default: None)
      -y VALUES, --values VALUES
                            Force 'count' or 'normalized' data parsing. Empty to auto-detect. (default:
                            None)
      -w, --cumm-levels     Activate if input table has already cummulative values on parent taxonomic
                            levels. (default: False)
      -s, --transpose       Transpose --input-table before parsing (if samples are listed on columns and
                            observations on rows) (default: False)
      -u [UNASSIGNED_HEADER ...], --unassigned-header [UNASSIGNED_HEADER ...]
                            Define one or more header names containing unsassinged/unclassified counts.
                            (default: None)
      -z REPLACE_ZEROS, --replace-zeros REPLACE_ZEROS
                            Treat zeros in the input table. INT (add 'smallest count' divided by INT to
                            every value), FLOAT (add FLOAT to every value). Default: 1000 (default:
                            1000)
      --obs-replace [OBS_REPLACE ...]
                            Replace values on observations labels/headers (supports regex). Example: '_'
                            ' ' will replace underscore with spaces, '^.+__' '' will remove the matching
                            regex. Several pairs of instructions are supported. (default: [])
      --sample-replace [SAMPLE_REPLACE ...]
                            Replace values on sample labels/headers (supports regex). Example: '_' ' '
                            will replace underscore with spaces, '^.+__' '' will remove the matching
                            regex. Several pairs of instructions are supported. (default: [])
      --min-frequency MIN_FREQUENCY
                            Define minimum number/percentage of samples containing an observation to
                            keep the observation [values between 0-1 for percentage, >1 specific
                            number]. (default: None)
      --max-frequency MAX_FREQUENCY
                            Define maximum number/percentage of samples containing an observation to
                            keep the observation [values between 0-1 for percentage, >1 specific
                            number]. (default: None)
      --min-count MIN_COUNT
                            Define minimum number/percentage of counts to keep an observation [values
                            between 0-1 for percentage, >1 specific number]. (default: None)
      --max-count MAX_COUNT
                            Define maximum number/percentage of counts to keep an observation [values
                            between 0-1 for percentage, >1 specific number]. (default: None)

    Samples options:
      -j TOP_OBS_BARS, --top-obs-bars TOP_OBS_BARS
                            Number of top abundant observations to show in the Samples panel, based on
                            the avg. percentage counts/sample. (default: 20)

    Heatmap and clustering options:
      -a {none,norm,log,clr}, --transformation {none,norm,log,clr}
                            Transformation of counts for Heatmap. none (counts), norm (percentage), log
                            (log10), clr (centre log ratio). (default: log)
      -e METADATA_COLS, --metadata-cols METADATA_COLS
                            Available metadata cols to be selected on the Heatmap panel. Higher values
                            will slow down the report navigation. (default: 3)
      --optimal-ordering    Activate optimal_ordering on scipy linkage method, takes longer for large
                            number of samples. (default: False)
      --show-zeros          Do not skip zeros on heatmap plot. File will be bigger and iteraction with
                            heatmap slower. By default, zeros will be omitted. (default: False)
      --linkage-methods [{single,complete,average,centroid,median,ward,weighted} ...]
      --linkage-metrics [{braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,kulczynski1,mahalanobis,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,yule} ...]
      --skip-dendrogram     Disable dendogram plots for clustering. (default: False)

    Correlation options:
      -x TOP_OBS_CORR, --top-obs-corr TOP_OBS_CORR
                            Number of top abundant observations to build the correlationn matrix, based
                            on the avg. percentage counts/sample. 0 for all (default: 50)


## Powered by

[<img src="https://static.bokeh.org/branding/logos/bokeh-logo.png" height="60">](https://bokeh.org)
[<img src="https://pandas.pydata.org/static/img/pandas.svg" height="40">](https://pandas.org)
[<img src="https://raw.githubusercontent.com/scipy/scipy/master/doc/source/_static/logo.svg" height="40">](https://scipy.org)
[<img src="http://scikit-bio.org/assets/logo.svg" height="40">](https://scikit-bio.org)
