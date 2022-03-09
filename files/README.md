# GRIMER References and other files

## Reference file format

1) File with a list (one per line) of taxonomic identifiers or taxonomic names

2) or formatted `.yml` file:

```yaml
"General Description":
  "Specific description":
    url: "www.website.com?id={}" 
    ids: [1,2,3]
```

The url can be a link to the entries listed on the id. Use the `{}` as a placeholder for the id. Example: `https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}`

The files should be provided in the main configuration file for grimer as follows:

```yaml
references:
  "Contaminants": "files/contaminants.yml"
  "Human-related": "files/human-related.yml" 
  "CUSTOM CONTAMINANTS": "file.txt"
  "LAB RELATED BACTERIA": "another_file.yml"
```

### contaminants.yml

Last update: 2022-03-09

Manually curated from diverse publications:

 | Organism group | Genus | Species | Reference |
 |----------------|-------|---------|-----------|
 | Bacteria | 6 | 0 | 1998 Tanner, M.A. et al. | 
 | Bacteria | 0 | 10 | 2002 Kulakov, L.A. et al. | 
 | Bacteria | 4 | 0 | 2003 Grahn, N. et al. | 
 | Bacteria | 16 | 0 | 2006 Barton, H.A. et al. | 
 | Bacteria | 11 | 1 | 2014 Laurence, M. et al.| 
 | Bacteria | 92 | 0 | 2014 Salter, S.J. et al. | 
 | Bacteria | 7 | 0 | 2015 Jervis-Bardy, J. et al. | 
 | Bacteria | 28 | 0 | 2015 Jousselin, E. et al. | 
 | Bacteria | 77 | 127 | 2016 Glassing, A. et al.| 
 | Bacteria | 23 | 0 | 2016 Lauder, A.P. et al. | 
 | Bacteria | 6 | 0 | 2016 Lazarevic, V. et al. | 
 | Bacteria | 62 | 0 | 2017 Salter, S.J. et al. | 
 | Bacteria | 0 | 122 | 2018 Kirstahler, P. et al. | 
 | Bacteria | 34 | 0 | 2018 Stinson, L.F. et al. | 
 | Bacteria | 18 | 0 | 2019 Stinson, L.F. et al. | 
 | Bacteria | 52 | 2 | 2019 Weyrich, L.S. et al. | 
 | Bacteria | 8 | 26 | 2019 de Goffau, M.C. et al. | 
 | Bacteria | 15 | 93 | 2020 Nejman D. et al. | 
 | Viruses | 0 | 1 | 2015 Kjartansdóttir, K.R. et al. | 
 | Viruses | 0 | 1 | 2015 Mukherjee, S. et al. | 
 | Viruses | 0 | 291 | 2019 Asplund, M. et al. |
 | Eukaryota | 0 | 3 | 2016 Czurda, S. et al. | 
 | Eukaryota | 0 | 1 | PRJNA168|
 | Total (unique) | 210 | 627 |  | 

### human-related.yml

Last update: 2022-03-09

Manually curated from from: Byrd, A., Belkaid, Y. & Segre, J. The human skin microbiome. Nat Rev Microbiol 16, 143–155 (2018). https://doi.org/10.1038/nrmicro.2017.157

```yaml
"Top organisms form the human skin microbiome":
  "Bacteria":
    url: "https://doi.org/10.1038/nrmicro.2017.157"
    ids: [257758, 225324, 169292, 161879, 146827, 43765, 38304, 38287, 38286, 29466, 29388, 28037, 1747, 1305, 1303, 1290, 1282, 1270]
  "Eukarya":
    url: "https://doi.org/10.1038/nrmicro.2017.157"
    ids: [2510778, 1047171, 379413, 119676, 117179, 76777, 76775, 76773, 44058, 41880, 36894, 34391, 31312, 5480, 5068, 3074, 2762]
  "Viruses":
    url: "https://doi.org/10.1038/nrmicro.2017.157"
    ids: [185639, 746832, 10566, 493803, 10279, 746830, 746831, 46771]
```

BacDive and eHOMD specific subsets. Dump date: 2022-03-09

```bash
scripts/bacdive_download.py
scripts/ehomd_download.py
```

## MGnify

The downloaded MGnify database file should be provided in the main configuration file for grimer as follows:

```yaml
external:
  mgnify: "files/mgnify5989.tsv"
```
### mgnify.tsv

MGnify dump date: 2022-03-09 (latest study accession MGYS00005989)

```bash
seq -f "MGYS%08g" 256 5989 | xargs -P 24 -I {} scripts/mgnify_download.py -i {} -v -g -o mgnify_dump_5989/ > mgnify_dump_5989.log 2>|1 |
scripts/mgnify_extract.py -f mgnify_dump_5989 -t 10 -o files/mgnify.tsv
```
