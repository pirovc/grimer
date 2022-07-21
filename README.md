# GRIMER

![GRIMER](grimer/img/logo.png)

GRIMER perform analysis of microbiome data and generates a portable and interactive dashboard integrating annotation, taxonomy and metadata.

## Examples

Online examples of reports generated with GRIMER: https://pirovc.github.io/grimer-reports/

## Installation

```bash
git clone https://github.com/pirovc/grimer.git
cd grimer
conda env create -f env.yaml # or mamba env create -f env.yaml
conda activate grimer # or source activate grimer
python setup.py install --record files.txt # Uninstall: xargs rm -rf < files.txt
grimer -h
```
***Soon GRIMER will be available as a package in BioConda.***

## Usage

### Tab-separated input table
```bash
grimer -i input_table.tsv
```

### BIOM file
```bash
grimer -i myfile.biom
```

### Tab-separated input table with taxonomic annotated observations (e.g. sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria...)
```bash
grimer -i input_table.tsv -f ";"
```

### Tab-separated input table with metadata
```bash
grimer -i input_table.tsv -m metadata.tsv
```

### With taxonomy integration (ncbi)
```bash
grimer -i input_table.tsv -m metadata.tsv -t ncbi #optional -b taxdump.tar.gz
```

### With configuration file to setup external tools, references and annotations
```bash
grimer -i input_table.tsv -m metadata.tsv -t ncbi -c config/default.yaml -d -g
```

### List all options 
```bash
grimer -h
```

### Analyzing any MGnify public study

```bash
./grimer-mgnify.py -i MGYS00006024 -o output_folder/
```

## Powered by

[<img src="https://static.bokeh.org/branding/logos/bokeh-logo.png" height="60">](https://bokeh.org)
[<img src="https://pandas.pydata.org/static/img/pandas.svg" height="40">](https://pandas.org)
[<img src="https://raw.githubusercontent.com/scipy/scipy/master/doc/source/_static/logo.svg" height="40">](https://scipy.org)
[<img src="http://scikit-bio.org/assets/logo.svg" height="40">](https://scikit-bio.org)
