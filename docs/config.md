# Configuration file

GRIMER uses a configuration file to set reference sources of annotation (e.g. contaminants), controls and external tools (decontam, mgnify). The configuration can be provided with the argument `-c/--config` and it should be in the [YAML](https://yaml.org/){ target="_blank" } format.

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

References can be provided as an external `.yml/.yaml` file in a specific format (see below) or as a text file with one taxonomic identifier or taxonomic name per line.

```yaml
"General Description":
  "Specific description":
    url: "www.website.com?id={}" 
    ids: [1,2,3]
```

A real example of saliva organisms extracted from BacDive (NCBI taxonomic ids):

```yaml
"Human-related bacterial isolates from BacDive":
  "Saliva":
    url: "https://bacdive.dsmz.de/search?search=taxid:{}"
    ids: [152331, 113107, 157688, 979627, 45634, 60133, 157687, 1624, 1583331, 1632, 249188]
```

Common contaminants compiled from the literature and human-related possible sources of contamination are available in the [GRIMER repository](https://github.com/pirovc/grimer/tree/main/files){ target="_blank" }. For more information, please refer to the [pre-print](https://doi.org/10.1101/2021.06.22.449360){ target="_blank" }. If the target study overlaps with some of those annotation (e.g. study of human skin), related entries can be easily removed from the provided files to not generate redundant annotations.

## controls

Several control groups can be provided to annotate samples. They can be provided as a file with one sample identifier per line:

```yaml
controls:
  "Controls": "controls.txt"
```

or directly from the metadata (`-m/--metadata-file`) as a field and value(s) information:

```yaml
controls:
  "Other Controls": 
    "sample_type": #  field
      - "blank"    #  value
      - "control"  #  value
```

Both methods can be combined into one configuration file.

## external

Set the configuration and functionality of external tools executed by GRIMER.

### mgnify

GRIMER uses a parsed MGnify database to annotate observations and link them to the respective MGnify repository, reporting most common biome occurrences. Instructions on how to re-generate the parsed database from MGnify can be found [here](https://github.com/pirovc/grimer/tree/main/files#mgnify){ target="_blank" }.

A [pre-parsed database](https://raw.githubusercontent.com/pirovc/grimer/main/files/mgnify5989.tsv){ target="_blank" } is available in the GRIMER repository (generated on 2022-03-09). To use it, please set the file in the configuration as follows and activate it with the `-g/--mgnify` when running GRIMER.

```yaml
external:
  mgnify: "files/mgnify5989.tsv"
```

### decontam

GRIMER can run [DECONTAM](https://benjjneb.github.io/decontam/){ target="_blank" } with `-d/--decontam`, but some configuration is necessary. It is possible to set the threshold (P* hyperparameter) and the method (frequency, prevalence, combined).

For the frequency/combined method, DNA frequencies for each sample have to be provided either in a `.tsv` file (sample identifier <tab> frequency) or as a metadata field. If none is provided, the sum of all counts in the input table is used for the frequency calculation.

For the prevalence/combined method, file(s) with a list of sample identifiers or a metadata field/value can be provided. If none is provided, all samples defined in the "controls" are considered for the prevalence calculation.

Below an example of how to set-up the configuration file for DECONTAM:

```yaml
external:
  decontam:
    threshold: 0.1 # P* hyperparameter threshold, values between 0 and 1
    method: "frequency" # Options: frequency, prevalence, combined
    frequency_file: "path/file1.txt"
    # frequency_metadata: "Field1"
    # prevalence_file: 
    #  - "path/file1.txt"
    #  - "path/file2.txt"
    prevalence_metadata: 
     "Field1":
      - "ValueA"
      - "ValueB"
      "Field2":
        - "ValueC"
```

## Using the configuration file

Example [UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom](https://microbiomedb.org/common/downloads/release-31/c66d2dc8473138e3a737ef2ad0b25f1e6e9c0f22/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom){ target="_blank" } file from [microbiomedb.org](https://microbiomedb.org){ target="_blank" }

config.yml (external .yml files are available in the [GRIMER repository](https://github.com/pirovc/grimer/tree/main/files){ target="_blank" })

```yml
references:
  "Contaminants": "files/contaminants.yml"
  "Human-related": "files/human-related.yml" 

external:
  mgnify: "files/mgnify5989.tsv"
  decontam:
    threshold: 0.1 # [0-1] P* hyperparameter
    method: "frequency" # frequency, prevalence, combined
```

Running GRIMER with DECONTAM and MGnify integration

```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom \
       --config config.yml \
       --decontam --mgnify \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species
```