# 04_pur_visualization — GenomeSyn

`run_genomesyn_*.sh` plot synteny with GenomeSyn 1.2.7:

- `run_genomesyn_denline_{5,10,100}kb.sh`, `run_genomesyn_denline_gdcolor.sh`: TE density lines
  (patched `GenomeSyn_densityline`), manual coords (`-cf1`), gene density (`-GD`), GD colouring.
- `run_genomesyn_TE.sh`, `run_genomesyn_TE_GD.sh`, `run_genomesyn_TE_GD_GENE2.sh`,
  `run_genomesyn_TE_GD_cf1.sh`: TE tracks with gene-density tracks (`-GD1`, `-GD2`).

Install GenomeSyn first, then set `GENOMESYN` in `../config.sh`.
The scripts activate the `genomesyn` environment, which provides the Perl runtime GenomeSyn needs.

## Environment

| Environment | Purpose |
|---|---|
| `genomesyn` | GenomeSyn 1.2.7 and the Perl runtime |

Environment file is `../envs/genomesyn.yml`.

Software versions: see the "Software versions" section of the top-level `README.md`.
