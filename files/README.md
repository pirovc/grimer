# GRIMER References and aux. files

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

Last update: 2021-04-01

 | Organism group | Genus | Species |
 |----------------|-------|---------|
 | Bacteria | 6 | 0 | 1998 Tanner, M.A. et al. |
 | Bacteria | 4 | 0 | 2003 Grahn, N. et al. |
 | Bacteria | 16 | 0 | 2006 Barton, H.A. et al. |
 | Bacteria | 11 | 1 | 2014 Laurence, M. et al. |
 | Bacteria | 92 | 0 | 2014 Salter, S.J. et al. |
 | Bacteria | 7 | 0 | 2015 Jervis-Bardy, J. et al. |
 | Bacteria | 28 | 0 | 2015 Jousselin, E. et al. | 
 | Bacteria | 23 | 0 | 2016 Lauder, A.P. et al. |
 | Bacteria | 6 | 0 | 2016 Lazarevic, V. et al. | 
 | Bacteria | 77 | 127 | 2016 Glassing, A. et al. |
 | Bacteria | 62 | 0 | 2017 Salter, S.J. et al. |
 | Bacteria | 0 | 122 | 2018 Kirstahler, P. et al. |
 | Bacteria | 8 | 26 | 2019 de Goffau, M.C. et al. | 
 | Bacteria | 52 | 2 | 2019 Weyrich, L.S. et al. |
 | Bacteria | 15 | 93 | 2020 Nejman D. et al. |
 | Viruses | 0 | 1 | 2015 Mukherjee, S. et al. |
 | Viruses | 0 | 1 | 2015 Kjartansdóttir, K.R. et al. |
 | Viruses | 0 | 301 | 2019 Asplund, M. et al. |
 | Total (unique) | 201 | 625 | |

### human-related.yml

BacDive and eHOMD dump date: 2021-04-13

```bash
scripts/bacdive_download.py
scripts/ehomd_download.py
```

## MGnify

The downloaded MGnify database file should be provided in the main configuration file for grimer as follows:

    external:
      mgnify: "files/mgnify.tsv"

## mgnify.tsv

MGnify dump date: 2021-04-08 (latest study accession MGYS00005724)

```bash
seq -f "MGYS%08g" 256 5724 | xargs -P 24 -I {} scripts/mgnify_download.py {} mgnify_dump_20210408/ > mgnify_dump_20210408.log 2>|1 |

scripts/mgnify_extract.py -f mgnify_dump_20210408 -t 10 -o files/mgnify.tsv
```