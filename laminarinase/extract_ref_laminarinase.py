#!/usr/bin/env python3
# coding: utf-8
import argparse
import csv
import brendapy.settings as settings
from brendapy import BrendaParser

"""
extract_ref_laminarinases.py

Extracts and annotates laminarinase accessions with EC classes and enzyme parameters (KM, pH optima, Topt) from BRENDA.
"""

def main():
    parser = argparse.ArgumentParser(
        description="Extract and annotate UniProt/PR accessions of laminarinases from BRENDA."
    )
    parser.add_argument(
        "--brenda-flatfile", "-b",
        default=settings.BRENDA_FILE,
        help="Path to BRENDA flat-file (defaults to bundled copy)"
    )
    parser.add_argument(
        "--ec-list", "-e",
        required=True,
        help="File listing EC numbers (one per line), e.g. laminarase_ec.txt"
    )
    parser.add_argument(
        "--out", "-o",
        default="laminarase_accessions_annotated.tsv",
        help="Output TSV file with accession annotations"
    )
    args = parser.parse_args()

    print(f"Loading BRENDA data from {args.brenda_flatfile} ...")
    brenda = BrendaParser(brenda_file=args.brenda_flatfile)

    ecs = [line.strip() for line in open(args.ec_list)
           if line.strip() and not line.startswith('#')]
    print("Found ECs:", ", ".join(ecs))

    metadata = {}  # accession -> {ecs:set, KM:list, pH:list, Topt:list}

    for ec in ecs:
        proteins = brenda.get_proteins(ec)
        print(f"Processing EC {ec}: {len(proteins)} entries")
        for prot in proteins.values():
            up = getattr(prot, 'uniprot', None)
            pr_ids = [d['data'] for d in prot.data.get('PR', [])]
            for acc in ([up] if up else []) + pr_ids:
                if not acc:
                    continue
                md = metadata.setdefault(acc, {'ecs': set(), 'KM': [], 'pH': [], 'Topt': []})
                md['ecs'].add(ec)
                for d in prot.data.get('KM', []):
                    md['KM'].append(str(d.get('value', '')))
                for d in prot.data.get('pH', []):
                    md['pH'].append(str(d.get('value', '')))
                for d in prot.data.get('Topt', []):
                    md['Topt'].append(str(d.get('value', '')))

    with open(args.out, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['accession', 'ec_numbers', 'KM_values', 'pH_optima', 'Topt'])
        for acc, md in sorted(metadata.items()):
            writer.writerow([
                acc,
                ';'.join(sorted(md['ecs'])),
                ';'.join(md['KM']),
                ';'.join(md['pH']),
                ';'.join(md['Topt']),
            ])
    print(f"Wrote {len(metadata)} annotated accessions to {args.out}")

if __name__ == '__main__':
    main()
