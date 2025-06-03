# Enzyme Finder

A reproducible pipeline for extracting enzyme data from BRENDA, fetching sequences from NCBI/UniProt, building proteome BLAST databases, and identifying homologs in a limited set of source organisms.

## Tools

- `extract-laminarinase`: extract annotated laminarinase entries
- `fetch-refs`: fetch protein sequences from NCBI/UniProt
- `get-taxid`: resolve taxonomic IDs
- `generate_proteome_db.sh`: download & construct BLAST databases
- `run_blast.sh`: run and filter BLAST hits

## Installation

```bash
pip install -e .
```

## 🧪 Example Data

The `examples/` folder contains both inputs and outputs for a full end-to-end run.

### Inputs
| File                             | Purpose                                 |
|----------------------------------|-----------------------------------------|
| `laminarinase_ec.txt`            | EC numbers to query BRENDA              |
| `species_taxids.tsv`             | Species and TaxIDs for proteome download |

### Outputs
| File                              | Description                                  |
|-----------------------------------|----------------------------------------------|
| `laminarase_accessions.tsv` | Extracted laminarinase annotations from BRENDA |
| `refs.fasta`                      | Protein sequences from NCBI/UniProt          |
| `all_proteomes.faa`              | Combined reference proteomes (for BLAST DB)  |
| `blast_raw.tsv`                  | Raw BLAST results                            |
| `blast_filtered.tsv`             | BLAST hits filtered by coverage, identity    |
| `blast_best_per_query.tsv`      | Best BLAST hit per query                     |


### Example Workflow

```bash
extract-laminarinase -e examples/laminarinase_ec.txt
fetch-refs -i laminarase_accessions.txt -e your@email
./bin/generate_proteome_db.sh examples/species_taxids.tsv output_dir
./bin/run_blast.sh refs.fasta output_dir/custom_proteomes output_dir

```
