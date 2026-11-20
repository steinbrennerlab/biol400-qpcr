# BIOL 400 qPCR runs

A website for browsing and downloading Bio-Rad CFX qPCR runs from BIOL 400. Every run gets
its own viewer page, where students select wells on a plate map and look at the amplification curves,
melt curves and melt peaks with Cq and Tm values. The index page lists all runs.

Run pages open with no wells selected. Wells with a Cq are outlined in blue, but their Cq and Tm
values only show once a student selects them. Students work out which wells are theirs from their
notes, select them, check the curves, and download the Cq/Tm CSV from the run page. Numbered
steps at the top of each run page say this. There is deliberately no "select all wells with a Cq"
button and no site-wide Cq export.

The site rebuilds itself (GitHub Actions → GitHub Pages) whenever anything in this repository
changes. Nobody needs to install anything to add a run.

## Adding a run

Make a new folder under `inputs/`, one folder per run, and name it starting with the date:

```
inputs/2026-11-17_lab3_AA/
    <run file>      the run data (see below)
    layout.csv      optional: sample names for each well
    run.json        optional: title, lab, section, machine, notes
```

**On github.com (no git needed):** open `inputs/`, choose **Add file → Upload files**, and drag
the whole run folder from your computer into the page. Then click **Commit changes**. The site updates
in about a minute. Progress is shown in the **Actions** tab.

### Run file: one of these

| File | Where it comes from |
|---|---|
| `.pcrd` | The data file CFX Manager / Maestro saves. Use it as is. It includes sample names if they were entered in CFX. |
| `.zpcr` | The raw file the instrument saves. Use it as is. |
| `.rdml` | CFX Maestro: *File → Export → Export RDML File*. Includes sample names. |
| Bio-Rad CSVs | CFX Maestro: *Export → Export All Data Sheets* (CSV). Put the CSVs in the run folder. |

If a folder has more than one (for example a `.pcrd` and the instrument's `.zpcr` for the same run),
the `.pcrd` is read and both are offered for download.

`.pcrd` files are encrypted. The build opens them with a password kept in the repository secret
`PCRD_PASSWORD` (Settings → Secrets and variables → Actions). The password is never stored in this
repository. If the secret is missing, `.pcrd`-only runs are listed as "pending".

### `layout.csv` (optional)

Sample names for each well. Only `Well` is required:

```
Well,Sample,Target,Content
A1,Leaf 1,ACTIN,Unkn
G1,water,PR1,NTC
```

`Content` uses Bio-Rad's well types (Unkn, NTC, Std, Pos Ctrl, Neg Ctrl, NRT). To get a blank
96-well template to fill in:
`python tools/qpcr_viewer.py inputs/<run>/<file>.zpcr --write-layout layout.csv`

### `run.json` (optional)

```json
{"title": "Lab 3: PR1 induction", "lab": "3", "section": "AA", "notes": "Plate sealed late"}
```

Without it, the folder name is used as the title and the date comes from the run file.
`lab`, `section` and the machine become filters on the index page.

### Machines

Each run shows which qPCR machine it ran on. The build reads the instrument serial from the run
file (or, if the file can't be opened, from its name) and looks it up under `"machines"` in
`site.json`:

```json
"machines": {"CC016428": "Imaizumi right machine", "BR006315": "Imaizumi left machine"}
```

To name a new machine, add its serial there. To override the machine for one run, put
`"machine": "..."` in its `run.json`.

## Previewing locally

```bash
python tools/build_site.py --serve     # builds _site/ and serves it at http://localhost:8000
```

To read `.pcrd` runs locally, set the password first (PowerShell: `$env:PCRD_PASSWORD = "..."`).

If a run can't be read, the build still finishes. That run is listed on the site with the reason,
and on GitHub it also appears as a warning on the Actions run. Use `--strict` to make the build
fail instead.

## How it's put together

| File | Role |
|---|---|
| `tools/qpcr_viewer.py` | Reads .pcrd / .zpcr / RDML / Bio-Rad CSV into one data structure and writes a standalone viewer page. It can also be run alone: `python tools/qpcr_viewer.py RUN.pcrd` |
| `tools/viewer_template.html` | The run viewer page (no external dependencies) |
| `tools/analysis.js` | Baseline, threshold, Cq, −dF/dT, Tm. **The only copy of the analysis.** It's inlined into every viewer page and loaded by the index for its per-run Cq counts, so the two always agree |
| `tools/build_site.py` | Turns `inputs/*` into `_site/`: index, viewer pages, run JSON, raw files and RFU CSVs |
| `tools/index_template.html` | The index page |
| `site.json` | Site title, intro text, and machine names by serial |
| `.github/workflows/pages.yml` | Builds and publishes the site on every push to `main` |

Only the Python standard library is needed (`openpyxl` only for `.xlsx` exports).
All analysis runs in the browser, so students can change the baseline and threshold and see the effect.

### One-time GitHub setup

In the repository: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
