# Data licences

The MIT licence in `LICENSE` covers the code, the documentation and the derived
tables in this repository. It does not cover the upstream data those tables are
derived from.

No upstream sequence, structure or annotation data is redistributed here.
`data/MANIFEST.tsv` records the source URL, the release and the SHA256 of each
input, and `make data` refetches it from the original providers, each under its
own licence:

| source | what is used | licence |
|---|---|---|
| UniProtKB | sequences, annotations, maturation features | CC BY 4.0 |
| AlphaFold DB | predicted structures | CC BY 4.0 |
| RCSB PDB | experimental structures | CC0 1.0 |
| Human Protein Atlas | expression, clusters | CC BY 4.0 |
| IEDB | B-cell assay records | as declared per record |
| PubMed / Europe PMC | retrieved text and bibliographic records | as declared per record |

Derived tables in `docs/results/` are redistributed under MIT; they contain
accession numbers and coordinates, not the upstream sequences.
