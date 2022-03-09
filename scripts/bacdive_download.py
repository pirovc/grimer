#!/usr/bin/env python3
import pandas as pd
from multitax import NcbiTx
import sys

## TODO
## filter infection?
## find names?

data = {("Host_Human-HostBodySite_Limb", "Limbs"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=39&filters%5B1%5D%5Bcat3%5D=&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodySite_Organ_Ear", "Ear"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=40&filters%5B1%5D%5Bcat3%5D=209&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodySite_Organ_Eye", "Eye"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=40&filters%5B1%5D%5Bcat3%5D=210&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodySite_Organ_Nose", "Nose"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=40&filters%5B1%5D%5Bcat3%5D=217&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodySite_Organ_SkinNailHair", "Skin/Nail/Hair"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=40&filters%5B1%5D%5Bcat3%5D=219&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodySite_Organ_OralCavityAndAirways", "Oral"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=5&filters%5B1%5D%5Bcat2%5D=41&filters%5B1%5D%5Bcat3%5D=&filters%5B1%5D%5Bcolor%5D=5&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1",
        ("Host_Human-HostBodyProduct_OralCavityAndAirways_Saliva", "Saliva"): "https://bacdive.dsmz.de/isolation-sources/csv?filter-domain=&filter-phylum=&filter-class=&filter-ordo=&filter-family=&filter-genus=&filters%5B0%5D%5Bcat1%5D=4&filters%5B0%5D%5Bcat2%5D=29&filters%5B0%5D%5Bcat3%5D=&filters%5B0%5D%5Bcolor%5D=4&filters%5B1%5D%5Bcat1%5D=6&filters%5B1%5D%5Bcat2%5D=47&filters%5B1%5D%5Bcat3%5D=276&filters%5B1%5D%5Bcolor%5D=6&iso_3_country=&polygon-strain-ids=&sort_by_order=ASC&sort_by=st.species&pfc=&csv=1"}

tax = NcbiTx(extended_names=True)

print('"Human-related bacterial isolates from BacDive":')

for (search, name), url in data.items():
    print('  "' + name + '":')
    print('    url: "https://bacdive.dsmz.de/search?search=taxid:{}"')
    parsed_ids = set()
    df = pd.read_table(url, sep=",", index_col=0).dropna(subset=["Species"])
    for species in df.Species.unique():
        taxids = tax.search_name(species, rank="species", exact=True)
        if not taxids:
            sys.stderr.write("Species name not found: " + species + "\n")
        elif len(taxids) > 1:
            sys.stderr.write("Species with ambiguous name: " + species + "\n")
        else:
            parsed_ids.add(taxids[0])
    print("    ids: [" + ", ".join(parsed_ids) + "]")
