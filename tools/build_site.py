#!/usr/bin/env python
"""Build the qPCR run browser site from the run folders in inputs/.

Each folder in inputs/ is one run:

    inputs/2026-11-17_lab3_AA/
        <run file>     .pcrd (needs PCRD_PASSWORD), .zpcr (raw instrument file), .rdml,
                       or Bio-Rad CSV/XLSX data sheets
        layout.csv     optional plate layout: Well, Sample, Target, Content (plate.csv also works)
        run.json       optional: {"title": "...", "date": "2026-11-17", "lab": "3", "section": "AA",
                                  "notes": "..."}

Output (default _site/):

    index.html                  table of runs, filters, combined Cq/Tm export
    analysis.js                 the Cq/Tm code (also inlined into every run page)
    runs/<id>.html              standalone viewer for that run (works offline once saved)
    runs/<id>.json              the run data the index uses for its combined export
    runs/<id>/...               raw input files and RFU tables as CSV

A run that can't be read is listed on the index with the reason, and the rest of the
site still builds (use --strict to fail instead). Only the Python standard library
is needed (plus openpyxl if a run uses .xlsx exports).

    python tools/build_site.py            # build into _site/
    python tools/build_site.py --serve    # build, then preview at http://localhost:8000
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import qpcr_viewer as qv  # noqa: E402

REPO = HERE.parent
INDEX_TEMPLATE = HERE / "index_template.html"
MARKER = ".qpcr-site"          # written into the output dir; only such a dir is ever wiped
LAYOUT_NAMES = ("layout.csv", "plate.csv", "layout.xlsx", "plate.xlsx")
RUN_SUFFIXES = (".pcrd", ".zpcr", ".rdml", ".csv", ".xlsx")
SINGLE_FILE_FORMATS = (".pcrd", ".zpcr", ".rdml")
IN_ACTIONS = os.environ.get("GITHUB_ACTIONS") == "true"


def warn(msg: str, title: str = "qPCR site"):
    if IN_ACTIONS:   # shows as an annotation on the workflow run page
        print(f"::warning title={title}::{msg.replace(chr(10), ' ')}")
    else:
        print(f"warning: {title}: {msg}", file=sys.stderr)


def slugify(name: str) -> str:
    s = re.sub(r"[^\w.-]+", "_", name).strip("._")
    return s or "run"


def read_run_meta(folder: Path) -> dict:
    p = folder / "run.json"
    if not p.exists():
        return {}
    try:
        meta = json.loads(p.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as e:
        raise qv.InputError(f"run.json is not valid JSON ({e})") from e
    if not isinstance(meta, dict):
        raise qv.InputError("run.json must be a JSON object like {\"title\": \"...\"}")
    return {k: str(v) for k, v in meta.items() if v is not None}


def find_layout(folder: Path) -> Path | None:
    for n in LAYOUT_NAMES:
        for p in folder.iterdir():
            if p.name.lower() == n:
                return p
    return None


def run_files(folder: Path, layout: Path | None) -> list[Path]:
    """Data files in the folder (and its subfolders, for Bio-Rad export folders)."""
    out = []
    for p in sorted(folder.rglob("*")):
        if not p.is_file() or p == layout or p.name.lower() == "run.json" or p.name.startswith("."):
            continue
        out.append(p)
    return out


def guess_date(meta: dict, run: dict | None, folder: Path) -> str:
    # run.json, then the folder name (local date, as people name runs), then the run file's
    # timestamp, which is UTC: an evening run in Seattle carries the next day's date.
    if meta.get("date"):
        return meta["date"]
    m = re.match(r"(\d{4}-\d{2}-\d{2})", folder.name)
    if m:
        return m.group(1)
    if run and run["meta"].get("run_date"):
        return run["meta"]["run_date"][:10]
    return ""


def write_rfu_csvs(run: dict, dest: Path, base: str) -> list[dict]:
    """Wide RFU tables (one row per cycle or temperature, one column per well) per channel."""
    wells = [qv.well_name(r, c) for r in range(run["plate"]["rows"]) for c in range(run["plate"]["cols"])]
    files = []
    for ch in run["channels"]:
        for kind, xs, xlabel, label in (("amp", run["cycles"], "Cycle", "Amplification RFU"),
                                        ("melt", run["melt_temps"], "Temperature", "Melt RFU")):
            curves = ch[kind]
            if not xs or not curves:
                continue
            cols = [w for w in wells if w in curves]
            suffix = "amplification" if kind == "amp" else "melt"
            name = f"{base}_{suffix}_{slugify(ch['name'])}.csv"
            with (dest / name).open("w", newline="", encoding="utf-8") as f:
                wr = csv.writer(f)
                wr.writerow([xlabel] + cols)
                for i, x in enumerate(xs):
                    wr.writerow([x] + ["" if curves[w][i] is None else curves[w][i] for w in cols])
            chan = f" ({ch['name']})" if len(run["channels"]) > 1 else ""
            files.append({"label": f"{label}{chan} CSV", "href": f"runs/{base}/{name}"})
    return files


def build_run(folder: Path, out: Path) -> dict:
    rid = slugify(folder.name)
    entry = {"id": rid, "folder": folder.name, "title": folder.name, "date": "", "lab": "",
             "section": "", "notes": "", "status": "ok", "message": "", "downloads": []}
    meta = read_run_meta(folder)
    entry.update({k: meta[k] for k in ("title", "lab", "section", "notes") if meta.get(k)})
    layout = find_layout(folder)
    files = run_files(folder, layout)
    readable = [p for p in files if p.suffix.lower() in RUN_SUFFIXES]

    if not readable:
        entry["date"] = guess_date(meta, None, folder)
        entry["status"] = "pending"
        entry["message"] = "No run file (.pcrd, .zpcr, .rdml or Bio-Rad CSV export) in this folder yet."
        return entry
    try:
        run, src = qv.load(readable)
    except qv.PcrdPasswordMissing:
        entry["date"] = guess_date(meta, None, folder)
        entry["status"] = "pending"
        entry["message"] = ("This run is a .pcrd file, and the build has no PCRD_PASSWORD to open it "
                            "(set it as a repository secret, or as an environment variable locally).")
        return entry
    if layout:
        qv.apply_layout(run, layout)
    title = entry["title"]
    run["meta"]["title"] = title
    run["meta"]["back_url"] = "../index.html"
    if entry["notes"]:
        run["meta"]["notes"].append(entry["notes"])
    entry["date"] = guess_date(meta, run, folder)

    runs_dir = out / "runs"
    dl_dir = runs_dir / rid
    dl_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / f"{rid}.html").write_text(qv.build_html(run, title), encoding="utf-8")
    # allow_nan=False: JSON.parse rejects NaN, so fail here rather than in the browser
    (runs_dir / f"{rid}.json").write_text(
        json.dumps(run, separators=(",", ":"), ensure_ascii=False, allow_nan=False), encoding="utf-8")

    # Offer every raw run file in the folder (e.g. both the .pcrd and the instrument's .zpcr),
    # the one that was read first.
    single = [p for p in readable if p.suffix.lower() in SINGLE_FILE_FORMATS]
    single.sort(key=lambda p: p != src)
    tables = [] if src.suffix.lower() in SINGLE_FILE_FORMATS else \
        [p for p in readable if p.suffix.lower() in (".csv", ".xlsx")]
    for p in single + tables + ([layout] if layout else []):
        shutil.copy2(p, dl_dir / p.name)
    downloads = [{"label": f"Run file (.{p.suffix.lstrip('.').lower()})", "href": f"runs/{rid}/{p.name}"} for p in single]
    downloads += [{"label": p.name, "href": f"runs/{rid}/{p.name}"} for p in tables]
    if layout:
        downloads.append({"label": "Plate layout CSV", "href": f"runs/{rid}/{layout.name}"})
    downloads += write_rfu_csvs(run, dl_dir, rid)
    downloads.append({"label": "Run data JSON", "href": f"runs/{rid}.json"})

    named = sum(1 for w in run["wells"].values() if w.get("sample"))
    entry.update(
        viewer=f"runs/{rid}.html", data=f"runs/{rid}.json", downloads=downloads,
        source_format=run["meta"].get("source_format", ""), source_file=src.name,
        channels=[c["name"] for c in run["channels"]],
        wells=run["plate"]["rows"] * run["plate"]["cols"], named_wells=named,
        has_layout=bool(layout or named),
    )
    return entry


def clean_out(out: Path):
    if out.exists():
        if not (out / MARKER).exists():
            raise SystemExit(f"refusing to overwrite {out}: it exists but wasn't made by this script "
                             f"(no {MARKER} file). Choose another --out or remove it yourself.")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    (out / MARKER).write_text("generated by tools/build_site.py; safe to delete\n", encoding="utf-8")
    (out / ".nojekyll").write_text("", encoding="utf-8")


def build(inputs: Path, out: Path, strict: bool = False) -> int:
    clean_out(out)
    site = {"title": "BIOL 400 qPCR runs", "intro": ""}
    site_file = REPO / "site.json"
    if site_file.exists():
        site.update(json.loads(site_file.read_text(encoding="utf-8-sig")))
    repo = os.environ.get("GITHUB_REPOSITORY")
    if repo:
        site["repo_url"] = f"https://github.com/{repo}"

    entries, errors = [], 0
    folders = sorted(p for p in inputs.iterdir() if p.is_dir() and not p.name.startswith((".", "_")))
    ids = set()
    for folder in folders:
        try:
            e = build_run(folder, out)
        except Exception as ex:   # one bad run must not take the whole site down
            errors += 1
            msg = str(ex) if isinstance(ex, qv.InputError) else f"{type(ex).__name__}: {ex}"
            warn(msg, title=f"inputs/{folder.name}")
            e = {"id": slugify(folder.name), "folder": folder.name, "title": folder.name,
                 "date": guess_date({}, None, folder), "lab": "", "section": "", "notes": "",
                 "status": "error", "message": f"Couldn't read this run: {msg}", "downloads": []}
            try:
                e.update({k: v for k, v in read_run_meta(folder).items() if k in ("title", "lab", "section", "date")})
            except qv.InputError:
                pass
        if e["id"] in ids:
            raise SystemExit(f"two run folders map to the same id '{e['id']}'; rename one of them")
        ids.add(e["id"])
        if e["status"] == "pending":
            warn(e["message"], title=f"inputs/{folder.name}")
        entries.append(e)
        print(f"  {e['status']:7}  {folder.name}" + (f"  ({e['message']})" if e["message"] else ""))

    entries.sort(key=lambda e: (e["date"], e["title"]), reverse=True)
    site["runs"] = entries
    site["built"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    shutil.copy2(HERE / "analysis.js", out / "analysis.js")
    payload = json.dumps(site, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    html = INDEX_TEMPLATE.read_text(encoding="utf-8")
    html = html.replace("__SITE_TITLE__", qv._html_escape(site["title"])).replace("__SITE_DATA__", payload)
    (out / "index.html").write_text(html, encoding="utf-8")

    ok = sum(e["status"] == "ok" for e in entries)
    print(f"built {out}: {ok} runs viewable, {len(entries) - ok} pending or with errors")
    return 1 if strict and errors else 0


def serve(out: Path, port: int):
    import functools
    import http.server
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(out))
    with http.server.ThreadingHTTPServer(("127.0.0.1", port), handler) as httpd:
        print(f"serving {out} at http://localhost:{port}/  (Ctrl+C to stop)")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", type=Path, default=REPO / "inputs", help="folder of run folders (default: inputs/)")
    ap.add_argument("--out", type=Path, default=REPO / "_site", help="output folder (default: _site/)")
    ap.add_argument("--strict", action="store_true", help="exit with an error if any run can't be read")
    ap.add_argument("--serve", action="store_true", help="after building, serve the site locally")
    ap.add_argument("--port", type=int, default=8000)
    a = ap.parse_args(argv)
    rc = build(a.inputs, a.out, a.strict)
    if a.serve:
        serve(a.out, a.port)
    return rc


if __name__ == "__main__":
    sys.exit(main())
