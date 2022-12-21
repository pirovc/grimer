# Importing files

GRIMER is independent of any quantification method and requires a contingency table with raw counts of observations/components for each samples/compositions in the study. Observations are usually, but not limited to, taxonomic entries (e.g. genus, species, strains), operational taxonomic units (OTUs), amplicon sequence variants (ASVs), metagenome-assembled genomes (MAGs) or sequence features.

GRIMER `--input-file` receives a tab-separated text file with a table of counts (Observation table, Count table, Contingency Tables, ...) or a [BIOM](https://biom-format.org/) file.

- Rows contain observations and columns contain samples (use `--transpose` if your file is reversed)
- First column and first row are used as headers
- Taxonomy integration: files can have either taxonomic identifiers (NCBI, e.g.: 562) or taxonomic names (NCBI, e.g.: Escherichia coli or GTDB, e.g.: s__Escherichia coli)

## biom file

GRIMER parses [BIOM](https://biom-format.org/) files and affiliated metadata (if available). Alternatively, an external metadata file can be provided with `-m/--metadata`.

Example [UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom](https://microbiomedb.org/common/downloads/release-31/c66d2dc8473138e3a737ef2ad0b25f1e6e9c0f22/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom) file from [microbiomedb.org](https://microbiomedb.org)

- Default report (no taxonomy)

```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom
```

- With integrated NCBI taxonomy (will translate names to taxids)

```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species
```
## tab-separated file (.tsv)

GRIMER parses .tsv files with single taxonomic identifier/names annotations or with multi-level (e.g.: lineage) taxonomic annotated observations.

### Multi-level annotations (e.g. Bacteria;Proteobacteria;Gammaproteobacteria...)

- Example [UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv](https://microbiomedb.org/common/downloads/release-31/c66d2dc8473138e3a737ef2ad0b25f1e6e9c0f22/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv) file from [microbiomedb.org](https://microbiomedb.org)


```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv \
       --level-separator ";"
```

- With metadata ([UgandaMaternalV3V4.16s_DADA2.sample_details.tsv](https://microbiomedb.org/common/downloads/release-31/c66d2dc8473138e3a737ef2ad0b25f1e6e9c0f22/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv))

```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv \
       --level-separator ";" \
       --metadata-file UgandaMaternalV3V4.16s_DADA2.sample_details.tsv
```

- With integrated NCBI taxonomy (will translate names to taxids)

```bash
grimer --input-file UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv \
       --level-separator ";" \
       --metadata-file UgandaMaternalV3V4.16s_DADA2.sample_details.tsv \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species
```

### Single level annotations (e.g. Neisseria animalis)

- Example [ERP108433_phylum_taxonomy_abundances_SSU_v4.1.tsv](https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00005180/pipelines/4.1/file/ERP108433_phylum_taxonomy_abundances_SSU_v4.1.tsv) from [MGnify](https://www.ebi.ac.uk/metagenomics), phylum level only

```bash
# Removing first column with kingdom
cut -f 2- ERP108433_phylum_taxonomy_abundances_SSU_v4.1.tsv > ERP108433_phylum_taxonomy_abundances_SSU_v4.1_parsed.tsv
# Set identifier for unassigned observations as "Unassigned" (many occurences, will be summed)
grimer --input-file ERP108433_phylum_taxonomy_abundances_SSU_v4.1_parsed.tsv \
       --unassigned-header "Unassigned"
```

- Re-generating taxonomic lineage from single annotations (in this case only superkingdom)

```bash
grimer --input-file ERP108433_phylum_taxonomy_abundances_SSU_v4.1_parsed.tsv \
       --unassigned-header "Unassigned" \
       --taxonomy ncbi \
       --ranks superkingdom phylum 
```

### Special cases
- --obs-replace --sample-replace --cumm-levels --transpose

## From commonly used tools/sources

### ganon

```bash
ganon table 
```

### MetaPhlAn
### QIIME2 feature table (.qza)

- Example [feature-table.qza](https://docs.qiime2.org/2022.8/data/tutorials/exporting/feature-table.qza) from [QIIME2 docs](https://docs.qiime2.org/2022.8/tutorials/exporting/#exporting-a-feature-table)

```bash
qiime tools export --input-path feature-table.qza --output-path exported-feature-table
grimer --input-file exported-feature-table/feature-table.biom
```
### MGnify
### phyloseq
### GTDB-tk
### CoverM
