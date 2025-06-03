# Enzyme Finder

A reproducible pipeline for extracting enzyme data from BRENDA, fetching sequences from NCBI/UniProt, building proteome BLAST databases, and identifying homologs.

## Tools

- `extract-laminarinase`: extract annotated laminarinase entries
- `fetch-refs`: fetch protein sequences from NCBI/UniProt
- `get-taxid`: resolve taxonomic IDs
- `generate_proteome_db.sh`: download & construct BLAST databases
- `run_blast.sh`: run and filter BLAST hits

## Installation

```bash
pip install -e .
