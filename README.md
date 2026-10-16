# MethaneWatch

Code and processed results for a test of methane plume attribution in the Permian Basin. The analysis uses Sentinel-5P TROPOMI observations from June to August 2024 and June to August 2026. It asks whether candidate enhancements are unusually close to gas flares after accounting for where oil and gas wells are located and where TROPOMI returned usable pixels.

[Results site](https://cpchoudhary2024.github.io/methanewatch/) | [Manuscript](paper/manuscript.md) | [Data sources](docs/DATA_SOURCES.md)

## Main result

| Summer | Overpasses | Candidates / days | Sign-inverted clusters | Flare enrichment: uniform control | Flare enrichment: well-density control | QA-pixel matched: uniform (95% day interval) | QA-pixel matched: wells (95% day interval) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2024 | 121 | 30 / 14 | 15 | 1.64 | 0.90 | 1.42 (0.93 to 1.83) | 0.86 (0.56 to 1.10) |
| 2026 | 119 | 10 / 5 | 16 | 1.97 | 1.08 | 1.69 (0.95 to 2.01) | 1.01 (0.56 to 1.24) |

Enrichment is the fraction of candidates within 5 km of a persistent flare divided by the corresponding control fraction. The uniform control samples the study box. The well-density control samples a kernel density built from 52,546 mapped wells. The matched analysis evaluates both controls on QA-passing pixels from each candidate's overpass and resamples whole candidate days. All four matched intervals include one. Proximity therefore does not isolate flares from surrounding infrastructure in these data.

Sign-inverted clusters are a detector diagnostic, not a calibrated false-discovery rate. Candidate locations and integrated mass enhancement rates are not facility identifications. Flares were mapped from 2019 to 2021, so the proximity test is spatial rather than a same-day activity test. The well layer comes from OpenStreetMap and may be incomplete. The next validation steps are listed in [`docs/DESIGN.md`](docs/DESIGN.md).

## Reproduce

Python 3.10 or newer is required. Downloading granules requires a free [Copernicus Data Space](https://dataspace.copernicus.eu/) account; copy `.env.example` to `.env` and add its credentials. The two cached seasons occupy about 14 GB.

```bash
make install
make test
make all
make opportunity
make serve
```

`make all` runs the 2026 Permian analysis and builds the manuscript and site. For the 2024 season, run:

```bash
make all START=2024-06-01 END=2024-09-01 OUT=outputs/permian-2024
```

`make opportunity` calculates the QA-pixel matched sensitivity for both seasons. `make serve` opens the 2026 results directory at `http://localhost:8000`. The other basin runs are archived in `outputs/`; their 2026 candidate counts were below the threshold for proximity tests. Input, output, and source hashes are recorded in each Permian season's `provenance.json`. Methods and fixed design constraints are in [DESIGN.md](docs/DESIGN.md).

## Repository

| Path | Contents |
|---|---|
| `src/methane/` | Detection, rates, spatial controls, mass balance, and site generation |
| `scripts/` | Figure, flare subset, and observation-matched analysis scripts |
| `tests/` | Offline unit tests |
| `outputs/` | Processed results and the generated site in `outputs/permian/` |
| `paper/` | Manuscript, figures, and submission assessment |

MIT license.
