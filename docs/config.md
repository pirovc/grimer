# Configuration file

GRIMER uses a configuration file to set reference sources of annotation (e.g. contaminants), controls and external tools (decontam, mgnify). The configuration can be provided with the argument `-c/--config` and it should be in the [YAML](https://yaml.org/) format.

A basic example of a configuration file:

```yaml
references:
  "Contaminants": "files/contaminants.yml"
  "Human-related": "files/human-related.yml" 

controls:
  "Negative Controls": "path/file1.tsv"
  "Positve Controls": 
    "Metadata_Field": 
      - "Metadata_Value1"
      - "Metadata_Value2"

external:
  mgnify: "files/mgnify5989.tsv"
  decontam:
    threshold: 0.1
    method: "frequency"
```

## references

References can be provided as an external `.yml/.yaml` file in a specific format (see below) or in a text file with one taxonomic identifier or taxonomic name per line.

```yaml
"General Description":
  "Specific description":
    url: "www.website.com?id={}" 
    ids: [1,2,3]
```

A real example of saliva organisms extracted from BacDive:

```yaml
"Human-related bacterial isolates from BacDive":
  "Saliva":
    url: "https://bacdive.dsmz.de/search?search=taxid:{}"
    ids: [152331, 113107, 157688, 979627, 45634, 60133, 157687, 1624, 1583331, 1632, 249188]
```

Common contaminants compiled from the literature and human-related possible sources of contamination are available in the [GRIMER repository](https://github.com/pirovc/grimer/tree/main/files). For more information, please refer to the [pre-print](https://doi.org/10.1101/2021.06.22.449360). If the target study overlaps with some of those annotation (e.g. study of human skin), related entries can be easily removed from those files to not generated redundant annotations.

## controls

Several control groups cann be provided. A file with one sample identifier per line

```yaml
controls:
  "Controls": "controls.txt"
```

or as a metadata field and value(s). 

```yaml
controls:
  "Other Controls": 
    "sample_type": #  field
      - "blank"    #  value
      - "control"  #  value
```

Both can be combined into one configuration file.

## external

Here it's possible to configure the function of external tools executed by GRIMER.

### mgnify

GRIMER can use a parsed MGnify database to annotate observations and link them to the respective MGnify entry, reporting most common biome occurrences. Instructions on how to re-generate the parsed database from MGnify can be found [here](https://github.com/pirovc/grimer/tree/main/files#mgnify).

A [pre-parsed database](https://raw.githubusercontent.com/pirovc/grimer/main/files/mgnify5989.tsv) is availble in the GRIMER repository (generated on 2022-03-09). To use it, please set the file in the configuration as follows and activate it with the `-g/--mgnify` when running GRIMER.

```yaml
external:
  mgnify: "files/mgnify5989.tsv"
```

### decontam

GRIMER can run [DECONTAM](https://benjjneb.github.io/decontam/) with `-d/--decontam`, but some configuration is necessary. It is possible to set the threshold (P* hyperparameter) and the method (frequency, prevalence, combined).

For the frequency/combined method, the DNA frequency for each sample has to be provided either in a `.tsv` file (sample identifier <tab> frequency) or a metadata field. If none is provided, the sum of all counts in the contigency table is used for the frequency calculation, but this is not recommended.

For the prevalence/combined method, file(s) with a list of sample identifiers can be provided or a field:value in the metadata. If none is provided, all samples defined in the "controls" are considered for the prevalence calculation.

Below an example of how to provide those values in the configuration file:

```yaml
external:
  decontam:
    threshold: 0.1 # [0-1] P* hyperparameter
    method: "frequency" # frequency, prevalence, combined
    # frequency_file: "path/file1.txt"
    # frequency_metadata: "Field1"
    # prevalence_file: 
    #  - "path/file1.txt"
    #  - "path/file2.txt"
    # prevalence_metadata: 
    #  "Field1":
    #    - "ValueA"
    #    - "ValueB"
    #  "Field2":
    #    - "ValueC"
```