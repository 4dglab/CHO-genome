# 07_run_EpiTAD-3D — 3D reconstruction

The scripts in this directory call EpiTAD-3D, which is not distributed with this repository.
Obtain it with:

```bash
git clone https://github.com/SaotomeMirai/EpiTAD-3D
```

| File | Purpose |
|---|---|
| `run_epiTAD.sh` | batch call of `bin/visualization/plot_singleLocus.py` for selected loci (10 kb resolution) |

Notes:

- `../envs/pastis.yml` recreates the environment from the EpiTAD-3D `environment.yml`
  (Python 3.9, cooler, numpy, pandas, matplotlib, scipy, scikit-learn, h5py, iced).
- Paths come from `../config.sh`: `EPITAD3D_HOME` (the clone), `EPITAD3D_MCOOL` (input matrix)
  and `EPITAD3D_OUT_DIR` (output root). Run `source config.sh` in the shell that starts the script.
- PASTIS-nb is bundled inside the EpiTAD-3D repository (`bin/third_party/`), so no extra install
  step is needed.

## Environment

| Environment | Purpose |
|---|---|
| `pastis` | cooler, numpy, pandas, matplotlib, scipy, scikit-learn, h5py, iced |

Environment file is `../envs/pastis.yml`.

Software versions: see the "Software versions" section of the top-level `README.md`.
