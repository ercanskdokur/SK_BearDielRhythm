# Example data — `BTR_example`

> **⚠️ SYNTHETIC DEMO DATA — NOT real animal locations or measurements.**
> These files were generated programmatically for one fictional individual
> (`BTR_example`) so that anyone can see the exact input formats the pipeline
> expects and run the first stages end-to-end. They contain no real GPS fixes,
> no real accelerometer readings, and no personal or sensitive information.

## What's here

| File | Rows | Purpose |
|------|-----:|---------|
| `Activity/BTR_example.txt` | 3600 (300 × 12) | 5-minute tri-axial accelerometer bins |
| `Position/BTR_example.txt` |  300 | ~hourly GPS fixes |
| `Bear_metadata.xlsx`       |    1 | one metadata row for `BTR_example` |

All three describe a single fictional bear tracked for ~12.5 days in May 2021,
inside the study bounding box. The filename stem (`BTR_example`) is the animal
ID the pipeline keys on, so the activity and position files **must share the
same stem**.

## How to use

Copy the files into the raw-input tree the pipeline reads (paths defined in
`00_setup.R`, relative to `PROJECT_ROOT/data`):

```
example_data/Activity/BTR_example.txt   ->  data/Activity/BTR_example.txt
example_data/Position/BTR_example.txt   ->  data/Position/BTR_example.txt
example_data/Bear_metadata.xlsx         ->  data/Bear_metadata.xlsx
```

Then set `PROJECT_ROOT` and run, inside the container:

```bash
PROJECT_ROOT=/path/to/project Rscript 01_data_preparation.R
PROJECT_ROOT=/path/to/project Rscript 02_response_variables.R
```

`01` → `02` run **offline** on this single bear (no GIS layers, no cluster) and
produce the cleaned activity/position objects and the hourly
intensity / duration response variables. Scripts from `03` onward require the
external covariate rasters (land cover, roads, dumps) and ERA5-Land temperature,
so they cannot be reproduced from the example alone — but the example still shows
exactly how every downstream script consumes the data.

## File formats (as parsed by `read_activity_file` / `read_position_file`)

Both readers are whitespace-delimited and **skip the first 3 header lines**.

**Activity** — 17 columns:
```
No CollarID UTC_Date UTC_Time LMT_Date LMT_Time Origin SCTS_Date SCTS_Time
Act_Mode DT ActivityX ActivityY ActivityZ Temp AnimalID GroupID
```
Dates are `DD/MM/YYYY`, times `HH:MM:SS`; `ActivityX/Y/Z` are integers in 0–255.
(The collars' LMT column equals UTC, hence LMT == UTC throughout the pipeline.)

**Position** — ≥14 whitespace fields per row; the parser reads field 5/6 as the
LMT date/time, field 7 as `Origin`, and detects latitude then longitude as the
first two **comma-decimal** numbers with ≥5 decimal places (e.g. `40,347677`).

**Metadata** — an Excel sheet with the columns
`Code, Name, Sex, Age (Year), Date, Mass (kg), Collar ID`
(`Date` as `DD/MM/YYYY`; `Code` must match the `.txt` filename stem).
