# GRIMER

![GRIMER](grimer/img/logo.png)

GRIMER automates analysis and reports an offline and interactive dashboard integrating annotation, taxonomy and metadata to analyse microbiome studies and detect contamination.

## Examples

Online examples of reports generated with GRIMER: https://pirovc.github.io/grimer-reports/

## Installation

```bash
git clone https://github.com/pirovc/grimer.git
cd grimer
conda env create -f env.yaml
conda activate grimer # source activate grimer
python setup.py install --record files.txt # Uninstall: xargs rm -rf < files.txt
grimer -h
```
***Soon GRIMER will be available as a package in BioConda.***

## Usage

### Basic
```bash
grimer -i input_table.tsv -m metadata.tsv -c config/default.yaml
```

### With taxonomy integration (ncbi)
```bash
grimer -i input_table.tsv -m metadata.tsv -c config/default.yaml -t ncbi #optional -b taxdump.tar.gz
```

### With DECONTAM and MGnify annotations
```bash
grimer -i input_table.tsv -m metadata.tsv -c config/default.yaml -d -g
```

### List all options 
```bash
grimer -h
```

## Powered by

[<img src="https://static.bokeh.org/branding/logos/bokeh-logo.png" height="150">](https://bokeh.org)
[<img src="https://pandas.pydata.org/static/img/pandas.svg" height="100">](https://pandas.org)
[<img src="https://raw.githubusercontent.com/scipy/scipy/master/doc/source/_static/logo.svg" height="100">](https://scipy.org)
[<img src="http://scikit-bio.org/assets/logo.svg" height="100">](https://scikit-bio.org)
