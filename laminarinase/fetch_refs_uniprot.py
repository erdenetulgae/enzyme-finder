#!/usr/bin/env python3
import time
import argparse
import os
from Bio import Entrez
import requests

def fetch_ncbi(acc, email):
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="protein", id=acc,
                               rettype="fasta", retmode="text")
        seq = handle.read()
        handle.close()
        if seq.startswith(">"):
            return seq
    except Exception:
        pass
    return None

def fetch_uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200 and r.text.startswith(">"):
            return r.text
    except Exception:
        pass
    return None

def main():
    p = argparse.ArgumentParser(
        description="Fetch FASTA for accession list, NCBI then UniProt fallback."
    )
    p.add_argument(
        "-i", "--accessions",
        default="laminarinase_accessions.txt",
        help="One accession per line (first token on each line)"
    )
    p.add_argument(
        "-o", "--out",
        default="data/refs.fasta",
        help="Output FASTA (path will be created if needed)"
    )
    p.add_argument(
        "-e", "--email",
        required=True,
        help="Your email for NCBI Entrez"
    )
    p.add_argument(
        "-d", "--delay",
        type=float,
        default=0.4,
        help="Seconds to wait between queries"
    )
    p.add_argument(
        "--skip-header",
        action="store_true",
        help="If set, skip the first line of the accessions file"
    )
    args = p.parse_args()

    # Ensure the output directory exists
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    Entrez.email = args.email

    # Read accessions, taking only the first whitespace-separated token on each line
    with open(args.accessions) as fh:
        accs = []
        for i, line in enumerate(fh):
            if i == 0 and args.skip_header:
                continue
            line = line.strip()
            if not line:
                continue
            # Split on whitespace and take first token
            acc = line.split()[0]
            accs.append(acc)

    success_count = 0
    with open(args.out, "w") as out_f:
        for acc in accs:
            print(f"Fetching {acc}...", end=" ", flush=True)
            seq = fetch_ncbi(acc, args.email) or fetch_uniprot(acc)
            if seq:
                out_f.write(seq)
                print("OK")
                success_count += 1
            else:
                print("FAILED")
            time.sleep(args.delay)

    print(f"\nDone: wrote {success_count} / {len(accs)} sequences to {args.out}")

if __name__ == "__main__":
    main()
