# GRIMER

<img src="https://raw.githubusercontent.com/pirovc/grimer/main/grimer/img/logo.png">

## About

GRIMER is a tool that performs automated analyses and generates a portable and interactive dashboard integrating annotation, taxonomy and metadata. It unifies several sources of evidence to help detect contamination. GRIMER is independent of quantification methods and directly analyses contingency tables to create an interactive and offline report. Reports can be created in seconds and are accessible for non-specialists, providing an intuitive set of charts to explore data distribution among observations and samples and its connections with external sources.

- More information about the method can be found in the [pre-print](https://doi.org/10.1101/2021.06.22.449360){ target="_blank" }
- The Source-code can be found in the [GitHub repository](https://github.com/pirovc/grimer){ target="_blank" }

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

- In-depth examples of input files can be found in [Importing files](importing)
- Complete examples of usage with real files can be found in [Examples](examples)


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

    grimer

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

    required arguments:
      -i INPUT_FILE, --input-file INPUT_FILE
                            Main input table with counts (Observation table, Count table, Contingency Tables, ...) or .biom file. By default rows contain observations and columns contain
                            samples (use --tranpose if your file is reversed). First column and first row are used as headers.

    main arguments:
      -m METADATA_FILE, --metadata-file METADATA_FILE
                            Input metadata file in simple tabular format with samples in rows and metadata fields in columns. QIIME 2 metadata format is also accepted, with an extra row to
                            define categorical and numerical fields. If not provided and --input-file is a .biom files, will attempt to get metadata from it.
      -t {ncbi,gtdb,silva,greengenes,ott}, --taxonomy {ncbi,gtdb,silva,greengenes,ott}
                            Define taxonomy to convert entry and annotate samples. Will automatically download and parse or files can be provided with --tax-files.
      -b [TAX_FILES ...], --tax-files [TAX_FILES ...]
                            Optional specific taxonomy files to use.
      -r [RANKS ...], --ranks [RANKS ...]
                            Taxonomic ranks to generate visualizations. Use 'default' to use entries from the table directly. Default: default
      -c CONFIG, --config CONFIG
                            Configuration file with definitions of references, controls and external tools.

    output arguments:
      -g, --mgnify          Plot MGnify chart
      -d, --decontam        Run and plot DECONTAM
      -l TITLE, --title TITLE
                            Title to display on the header of the report.
      -p [{overview,samples,heatmap,correlation} ...], --output-plots [{overview,samples,heatmap,correlation} ...]
                            Plots to generate. Default: overview,samples,heatmap,correlation
      -o OUTPUT_HTML, --output-html OUTPUT_HTML
                            File to output report. Default: output.html
      --full-offline        Embed javascript library in the output file. File will be around 1.5MB bigger but also work without internet connection. That way your report will live forever.

    general data options:
      -f LEVEL_SEPARATOR, --level-separator LEVEL_SEPARATOR
                            If provided, consider --input-table to be a hierarchical multi-level table where the observations headers are separated by the indicated separator characther
                            (usually ';' or '|')
      -y VALUES, --values VALUES
                            Force 'count' or 'normalized' data parsing. Empty to auto-detect.
      -w, --cumm-levels     Activate if input table has already cummulative values among levels.
      -s, --transpose       Transpose --input-table (if samples are listed on columns and observations on rows)
      -u [UNASSIGNED_HEADER ...], --unassigned-header [UNASSIGNED_HEADER ...]
                            Define one or more header names containing unsassinged/unclassified counts.
      --obs-replace [OBS_REPLACE ...]
                            Replace values on table observations labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching
                            regex.
      --sample-replace [SAMPLE_REPLACE ...]
                            Replace values on table sample labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.
      -z REPLACE_ZEROS, --replace-zeros REPLACE_ZEROS
                            INT (add 'smallest count'/INT to every raw count), FLOAT (add FLOAT to every raw count). Default: 1000
      --min-frequency MIN_FREQUENCY
                            Define minimum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].
      --max-frequency MAX_FREQUENCY
                            Define maximum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].
      --min-count MIN_COUNT
                            Define minimum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].
      --max-count MAX_COUNT
                            Define maximum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].

    Samples options:
      -j TOP_OBS_BARS, --top-obs-bars TOP_OBS_BARS
                            Top abundant observations to show in the bars.

    Heatmap and clustering options:
      -a TRANSFORMATION, --transformation TRANSFORMATION
                            none (counts), norm (percentage), log (log10), clr (centre log ratio). Default: log
      -e METADATA_COLS, --metadata-cols METADATA_COLS
                            How many metadata cols to show on the heatmap. Higher values makes plot slower to navigate.
      --optimal-ordering    Activate optimal_ordering on linkage, takes longer for large number of samples.
      --show-zeros          Do not skip zeros on heatmap. File will be bigger and iteraction with heatmap slower.
      --linkage-methods [{single,complete,average,centroid,median,ward,weighted} ...]
      --linkage-metrics [{braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,wminkowski,yule} ...]
      --skip-dendrogram     Disable dendogram. Will create smaller files.

    Correlation options:
      -x TOP_OBS_CORR, --top-obs-corr TOP_OBS_CORR
                            Top abundant observations to build the correlationn matrix, based on the avg. percentage counts/sample. 0 for all

## Powered by

[<img src="https://static.bokeh.org/branding/logos/bokeh-logo.png" height="60">](https://bokeh.org)
[<img src="https://pandas.pydata.org/static/img/pandas.svg" height="40">](https://pandas.org)
[<img src="https://raw.githubusercontent.com/scipy/scipy/master/doc/source/_static/logo.svg" height="40">](https://scipy.org)
[<img src="http://scikit-bio.org/assets/logo.svg" height="40">](https://scikit-bio.org)
