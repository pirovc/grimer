#!/usr/bin/env python3
import pandas as pd
import sys
import urllib.request
import re


def get_taxid(url):
    try:
        sys.stderr.write(url + "\n")
        assembly_stats = url + "/" + url.split("/")[-1] + "_assembly_stats.txt"
        filedata = urllib.request.urlopen(assembly_stats).read().decode()
        x = re.search("# Taxid:[\s0-9]*\\r\\n", filedata)
        if x:
            return re.findall("\d+", x.group())[0]
        else:
            return None
    except:
        return None

# Can be Oral, Nasal or both ("Nasal,Oral")
habitats = ["Oral", "Nasal"]
data = "http://www.ehomd.org/ftp/genomes/PROKKA/current/SEQID_info.csv"

df = pd.read_table(data, sep=",", usecols=["Habitat", "Sequence_Source"])
df = df[df["Habitat"].isin(habitats + ["Nasal,Oral"])].drop_duplicates()
df["taxid"] = df["Sequence_Source"].map(get_taxid)

print('"Human Oral Microbiome Database (eHOMD)":')
for h in habitats:
    print('  "' + h + '":')
    parsed_ids = set(df.taxid[df.Habitat.str.contains(h)])
    print('    url: "http://www.ehomd.org/?name=HOMD"')
    print("    ids: [" + ", ".join(parsed_ids) + "]")

sys.stderr.write("Could not retrieve taxid for: " + "\n".join(df[df.taxid.isna()]["Sequence_Source"].to_list()) + "\n")
