# Lexibank Forms with Alignments

Reduced extract of the consolidated Lexibank cross-linguistic wordlist
database. The original CSV (~100 MB, 21 columns) has been reduced to a
6-column TSV keeping only the fields needed for alignment experiments.

## Columns

| Column | Description |
|--------|-------------|
| Dataset | Source dataset identifier |
| Language_ID | Language identifier |
| Concepticon_Gloss | Standardized semantic concept |
| Segments | Space-separated phonetic segments |
| Cognacy | Cognate class identifier |
| Alignment | Gold-standard alignment (space-separated, with gaps) |

## Provenance

Aggregated from 58 datasets via Lexibank / CLDF pipelines. Original data
compiled by the respective dataset authors; see each dataset's metadata
for licensing details.
