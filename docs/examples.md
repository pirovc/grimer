![GRIMER](https://raw.githubusercontent.com/pirovc/grimer/main/grimer/img/logo.png)

Examples of reports generated with [GRIMER](https://github.com/pirovc/grimer)

---

### Data analysis from Leiby et al. "Lack of detection of a human placenta microbiome in samples from preterm and term deliveries"

***original publication: [10.1186/s40168-018-0575-4](https://doi.org/10.1186/s40168-018-0575-4){ target="_blank" }***

**[GRIMER report MGS](https://pirovc.github.io/grimer-reports/placenta/placenta_mgs.html){ target="_blank" }**

**[GRIMER report AMPLICON](https://pirovc.github.io/grimer-reports/placenta/placenta_amplicon.html){ target="_blank" }**

<details>
<summary>commands used to create report</summary>

```bash
# Download files (table, metadata and config)
wget https://raw.githubusercontent.com/pirovc/grimer-reports/main/placenta/placenta_files.tar.gz
tar xf placenta_files.tar.gz

# Run GRIMER
# AMPLICON
grimer --config placenta_amplicon_config.yaml \
       --input-file placenta_amplicon_table.tsv \
       --metadata-file placenta_metadata.tsv \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species \
       --level-separator ";" \
       --obs-replace "^.+__" "" "_" " " \
       --unassigned-header "Unassigned"  \
       --decontam --mgnify --transpose \
       --title "Placenta study AMPLICON - Leiby, J.S. et al 2018" \
       --output-html placenta_amplicon.html

# MGS
grimer --config placenta_mgs_config.yaml \
       --input-file placenta_mgs_table.tsv \
       --metadata-file placenta_metadata.tsv \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species \
       --level-separator "|" \
       --unassigned-header "unassigned"  \
       --decontam --mgnify \
       --title "Placenta study MGS - Leiby, J.S. et al 2018" \
       --output-html placenta_mgs.html
```
</details>

---

### KatharoSeq analysis from Minich et al. "KatharoSeq Enables High-Throughput Microbiome Analysis from Low-Biomass Samples"

***original publication: [10.1128/mSystems.00218-17](https://doi.org/10.1128/mSystems.00218-17){ target="_blank" }***

**[GRIMER report](https://pirovc.github.io/grimer-reports/katharoseq/katharoseq.html){ target="_blank" }**

<details>
<summary>commands used to create report</summary>

```bash
# Download files (table, metadata and config)
wget https://raw.githubusercontent.com/pirovc/grimer-reports/main/katharoseq/katharoseq_files.tar.gz
tar xf katharoseq_files.tar.gz

# Run GRIMER
grimer --config katharoseq_config.yaml \
       --input-file katharoseq_table.tsv \
       --metadata-file katharoseq_metadata.tsv \
       --transformation clr \
       --obs-replace "^.+__" "" "_" " " \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species \
       --level-separator ";" \
       --decontam --mgnify \
       --title "KatharoSeq - Minich et al. 2018" \
       --output-html katharoseq.html
```

</details>

---

### Preterm Infant Resistome downloaded from [MicrobiomeDB](https://microbiomedb.org/mbio/app/record/dataset/DS_82fe0308e2){ target="_blank" }

***original publication: [10.1038/nmicrobiol.2016.24](https://doi.org/10.1038/nmicrobiol.2016.24){ target="_blank" }***

**[GRIMER report](https://pirovc.github.io/grimer-reports/microbiomedb/ResistomeAmplicon.html){ target="_blank" }**

<details>
<summary>commands used to create report</summary>

```bash
# Download files (table, metadata and config) - Original source: https://microbiomedb.org/common/downloads/release-22/82fe0308e2032de2041694df6592ba542ea84b86/ResistomeAmplicon.16s_DADA2.taxon_abundance.biom
wget https://raw.githubusercontent.com/pirovc/grimer-reports/main/microbiomedb/microbiomedb_files.tar.gz
tar xf microbiomedb_files.tar.gz

# Run GRIMER
grimer --config ResistomeAmplicon.16s_DADA2_config.yaml \
       --input-file ResistomeAmplicon.16s_DADA2.taxon_abundance.biom \
       --taxonomy ncbi \
       --ranks superkingdom phylum class order family genus species \
       --decontam --mgnify \
       --title "MicrobiomeDB Preterm Infant Resistome (V4)" \
       --output-html ResistomeAmplicon.html
```

</details>

---

### Antibiotic induced changes in the microbiota disrupt redox dynamics in the gut downloaded from [MGnify](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005180){ target="_blank" }

***original publication [10.7554/elife.35987](https://doi.org/10.7554/elife.35987){ target="_blank" }***

**[GRIMER report](https://pirovc.github.io/grimer-reports/mgnify/MGYS00005180.html){ target="_blank" }**

<details>
<summary>commands used to create report</summary>

```bash
# Script to download files and generate GRIMER report from any MGnify study accession
# Requires "jsonapi-client>=0.9.7" (conda install "jsonapi-client>=0.9.7")
./grimer-mgnify.py -i MGYS00005180 -o MGYS00005180 -g "--decontam --mgnify" 

# Or directly from files
wget https://raw.githubusercontent.com/pirovc/grimer-reports/main/mgnify/mgnify_files.tar.gz
tar xf mgnify_files.tar.gz
# Run GRIMER
grimer --config MGYS00005180_config.yaml \
       --input-file MGYS00005180_ERP108433_taxonomy_abundances_SSU_v4.1.tsv \
       --metadata-file MGYS00005180_metadata.tsv \
       --obs-replace "^.+__" "" "_" " " \
       --taxonomy ncbi \
       --ranks superkingdom kingdom phylum class order family genus species \
       --level-separator ";" \
       --decontam --mgnify \
       --title "MGnify study accession MGYS00005180" \
       --output-html MGYS00005180.html
```

</details>

---