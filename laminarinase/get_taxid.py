#!/usr/bin/env python3
import re
from Bio import Entrez

Entrez.email = "erka.enkhbold23@imperial.ac.uk"

# regex to strip off common strain/collection identifiers + anything after
STRAIN_RE = re.compile(
    r"\b(?:ATCC|DSMZ?|NBRC|PCC|PAO|MG1655|KT2440|CN-32|UCC2003|BP-1|ES114|81-176|13882|40230)\b.*",
    flags=re.IGNORECASE,
)

def clean_species_name(full_name):
    """Strip off strain/collection codes so we get just 'Genus species'."""
    # remove e.g. 'Bacillus megaterium DSM319' → 'Bacillus megaterium'
    no_strain = STRAIN_RE.sub("", full_name)
    parts = no_strain.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return full_name

def fetch_taxid(species):
    """Try full name first, then fallback to cleaned genus+species."""
    # 1) full‐phrase search
    term = f'"{species}"[Scientific Name]'
    rec = Entrez.read(Entrez.esearch(db="taxonomy", term=term))
    if rec["IdList"]:
        return rec["IdList"][0]
    
    # 2) fallback to genus+species only
    base = clean_species_name(species)
    if base != species:
        rec2 = Entrez.read(Entrez.esearch(db="taxonomy", term=f'"{base}"[Scientific Name]'))
        if rec2["IdList"]:
            # if multiple, report them for manual curation
            if len(rec2["IdList"]) > 1:
                print(f"⚠ Multiple matches for '{base}': {rec2['IdList']}")
            return rec2["IdList"][0]
    
    return "NOT_FOUND"

def main():
    with open("data/species.txt") as inf, open("data/species_taxids.tsv","w") as out:
        out.write("species\ttaxid\n")
        for line in inf:
            sp = line.strip()
            taxid = fetch_taxid(sp)
            out.write(f"{sp}\t{taxid}\n")
    print("Done → species_taxids.tsv")

main()

