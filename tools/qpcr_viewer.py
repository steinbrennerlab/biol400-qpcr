#!/usr/bin/env python
"""Build a standalone HTML viewer for Bio-Rad CFX qPCR runs.

Students open the resulting .html file in any browser (no internet, no install)
and browse amplification and melt curves for any group of wells they select.

Supported inputs
----------------
* ``.zpcr``  raw run file written by a CFX96 / CFX Opus instrument (read directly).
* ``.rdml``  RDML export from CFX Maestro / CFX Manager (File > Export > RDML).
* Bio-Rad CSV / XLSX data-sheet exports ("Export All Data Sheets"), given as a
  folder or as individual files: Quantification Amplification Results,
  Melt Curve RFU Results, and optionally Quantification Cq Results /
  Quantification Summary / Melt Curve Peak Results for sample names.

``.pcrd`` files are password-protected archives that only Bio-Rad software can
open. Open the .pcrd in CFX Maestro and export RDML (or the data sheets), then
run this script on the export.

Usage
-----
    python qpcr_viewer.py RUN.zpcr                      # -> RUN.html
    python qpcr_viewer.py RUN.zpcr --layout plate.csv   # add sample names
    python qpcr_viewer.py RUN.zpcr --write-layout plate.csv   # blank layout to fill in
    python qpcr_viewer.py export_folder/ -o lab3.html --title "Lab 3 qPCR"

The layout CSV has columns Well, Sample, Target, Content (Content is Bio-Rad's
well type: Unkn, NTC, Std, Pos Ctrl, Neg Ctrl, NRT). Only Well is required.

Only the Python standard library is needed (openpyxl only for .xlsx exports).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import re
import struct
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "viewer_template.html"
ANALYSIS_JS = HERE / "analysis.js"   # Cq/Tm code shared with the site index; inlined into each page

# CFX96 optical channels, in the order they are stored in raw plate reads.
CFX96_CHANNELS = ["FAM/SYBR", "HEX/VIC", "Texas Red/ROX", "Cy5", "Quasar 705", "FRET"]
ROWS = "ABCDEFGHIJKLMNOP"


class InputError(Exception):
    pass


# --------------------------------------------------------------------------
# Common data model
# --------------------------------------------------------------------------
def well_name(r: int, c: int) -> str:
    return f"{ROWS[r]}{c + 1}"


def norm_well(s: str) -> str | None:
    m = re.fullmatch(r"\s*([A-Pa-p])\s*0*(\d{1,2})\s*", str(s))
    return f"{m.group(1).upper()}{int(m.group(2))}" if m else None


def new_run(rows=8, cols=12):
    return {
        "meta": {"notes": []},
        "plate": {"rows": rows, "cols": cols},
        "cycles": [],
        "melt_temps": [],
        "channels": [],          # [{name, amp:{well:[..]}, melt:{well:[..]}}]
        "wells": {},             # well -> {sample, target, content}
        "vendor": {},            # channel -> well -> {cq, tm}
    }


def get_channel(run, name):
    for ch in run["channels"]:
        if ch["name"] == name:
            return ch
    ch = {"name": name, "amp": {}, "melt": {}}
    run["channels"].append(ch)
    return ch


# --------------------------------------------------------------------------
# .zpcr  (raw instrument run file: an ordinary zip of binary plate reads)
# --------------------------------------------------------------------------
DATE_RE = re.compile(rb"[A-Z][a-z]{2}, \d{1,2} [A-Z][a-z]{2} \d{4} \d\d:\d\d:\d\d GMT")


def parse_plateread(b: bytes):
    """Decode one ReadNNNNN.Plateread record.

    Layout (reverse-engineered from CFX96 Touch files, firmware DSP 1.123):
      ... null-terminated firmware string ...
      7 big-endian uint32: 1, 3, read_no, 0, protocol_step (0-based), repeat_no, 0
      1 pad byte, 5 big-endian float32 (block temp, ..., lid temp)
      ... date string 'Tue, 17 Nov 2026 12:48:26 GMT' ...
      5 bytes, uint32 LE n_floats, then n_floats LE float32 =
      6 channels x (96 wells + 12 reference positions) x (mean, sd, min, max)
    """
    m = DATE_RE.search(b)
    if not m:
        raise InputError("plate read without timestamp; unrecognised .zpcr layout")
    p = m.end() + 5
    (n_floats,) = struct.unpack("<I", b[p:p + 4])
    p += 4
    vals = struct.unpack(f"<{n_floats}f", b[p:p + 4 * n_floats])
    recs = n_floats // 4
    per_ch = recs // 6
    if per_ch == 108:
        rows, cols = 8, 12
    else:
        raise InputError(f"unsupported plate layout ({per_ch} positions per channel); "
                         "only 96-well CFX runs are supported for raw .zpcr files")
    nwell = rows * cols
    chans = []
    for c in range(6):
        base = c * per_ch * 4
        chans.append([vals[base + 4 * w] for w in range(nwell)])

    hdr = {"date": m.group().decode()}
    fw = b.find(b"Firmware versions")
    if fw >= 0:
        null = b.find(b"\x00", fw)
        try:
            ints = struct.unpack(">7I", b[null + 2:null + 30])
            temps = struct.unpack(">5f", b[null + 31:null + 51])
            if ints[0] == 1 and ints[3] == 0 and ints[4] < 100 and 0 < ints[5] < 10000:
                hdr.update(read_no=ints[2], step=ints[4], repeat=ints[5], block_temp=temps[0])
        except struct.error:
            pass
    return hdr, chans, rows, cols


def parse_protocol(text: str):
    """Return list of TEMP steps: [{temp, hold, inc, read}] (0-based index = step id)."""
    steps = []
    for tok in text.split(";"):
        tok = tok.strip()
        if tok.startswith("TEMP "):
            t, h = tok[5:].split(",")[:2]
            steps.append({"temp": float(t), "hold": float(h), "inc": 0.0, "read": False})
        elif tok.startswith("INC ") and steps:
            steps[-1]["inc"] = float(tok[4:])
        elif tok.startswith("PLATEREAD") and steps:
            steps[-1]["read"] = True
    return steps


def describe_protocol(text: str) -> str:
    out, steps = [], 0
    for tok in text.split(";"):
        tok = tok.strip()
        if tok.startswith("TEMP "):
            steps += 1
            t, h = tok[5:].split(",")[:2]
            out.append(f"{float(t):g}°C {int(float(h))}s")
        elif tok.startswith("PLATEREAD"):
            out[-1] += " +read"
        elif tok.startswith("INC "):
            out[-1] += f" (+{tok[4:]}°C/step)"
        elif tok.startswith("GOTO "):
            n, r = tok[5:].split(",")
            out.append(f"repeat ×{int(r) + 1}")
    return " → ".join(out)


def read_zpcr(path: Path):
    z = zipfile.ZipFile(path)
    names = sorted(n for n in z.namelist() if n.lower().endswith(".plateread"))
    if not names:
        raise InputError(f"{path.name}: no plate reads found (run may not have started)")
    proto_txt = z.read("ProtocolRunDefinition.txt").decode("latin-1") \
        if "ProtocolRunDefinition.txt" in z.namelist() else ""
    steps = parse_protocol(proto_txt)

    info = {}
    if "RunInfo.xml" in z.namelist():
        root = ET.fromstring(z.read("RunInfo.xml"))
        for kv in root.iter("KeyValuePairs"):
            k, v = kv.findtext("Key"), kv.findtext("Value")
            if k:
                info[k] = (v or "").strip().strip('"')

    reads = [parse_plateread(z.read(n)) for n in names]
    rows, cols = reads[0][2], reads[0][3]
    run = new_run(rows, cols)

    # Group reads by protocol step. A step with an INC is a melt; the other read
    # step is the amplification (quantification) step.
    groups: dict[int, list] = {}
    if all("step" in h for h, *_ in reads):
        for h, ch, *_ in reads:
            groups.setdefault(h["step"], []).append((h, ch))
    else:  # fallback: assume reads are in order, split where block temp starts climbing
        raise InputError("could not decode plate-read headers in this .zpcr")

    amp_step = melt_step = None
    for s in groups:
        is_melt = s < len(steps) and steps[s]["inc"] > 0
        if is_melt and melt_step is None:
            melt_step = s
        elif not is_melt and amp_step is None:
            amp_step = s

    active = [c for c in range(6) if any(any(v != 0 for v in ch[c]) for _, ch in
                                         (groups.get(amp_step) or groups.get(melt_step)))]
    if amp_step is not None:
        g = sorted(groups[amp_step], key=lambda x: x[0]["repeat"])
        run["cycles"] = [h["repeat"] for h, _ in g]
        for c in active:
            chn = get_channel(run, CFX96_CHANNELS[c])
            for w in range(rows * cols):
                chn["amp"][well_name(w // cols, w % cols)] = [round(ch[c][w], 2) for _, ch in g]
    if melt_step is not None:
        g = sorted(groups[melt_step], key=lambda x: x[0]["repeat"])
        st = steps[melt_step]
        run["melt_temps"] = [round(st["temp"] + st["inc"] * (h["repeat"] - 1), 2) for h, _ in g]
        for c in active:
            chn = get_channel(run, CFX96_CHANNELS[c])
            for w in range(rows * cols):
                chn["melt"][well_name(w // cols, w % cols)] = [round(ch[c][w], 2) for _, ch in g]

    first_date = reads[0][0].get("date", "")
    try:
        run_date = dt.datetime.strptime(first_date, "%a, %d %b %Y %H:%M:%S GMT").strftime("%Y-%m-%d %H:%M UTC")
    except ValueError:
        run_date = first_date
    run["meta"].update(
        source_format="Bio-Rad raw run file (.zpcr)",
        run_date=run_date,
        instrument=f"Bio-Rad CFX ({info.get('BlockDescription', '96')} block, serial {info.get('BaseSerialNumber', '?')})",
        protocol=describe_protocol(proto_txt),
        rfu_note="Raw fluorescence as read by the instrument (no Bio-Rad dye calibration or "
                 "background correction applied).",
    )
    if "ended" not in z.namelist():
        run["meta"]["notes"].append("This run file has no 'ended' marker: the run may have been "
                                    "stopped early.")
    return run


# --------------------------------------------------------------------------
# RDML export
# --------------------------------------------------------------------------
SAMPLE_TYPES = {"unkn": "Unkn", "ntc": "NTC", "std": "Std", "pos": "Pos Ctrl",
                "neg": "Neg Ctrl", "nrt": "NRT", "ntp": "NTP", "opt": "Opt"}


def read_rdml(path: Path):
    raw = path.read_bytes()
    if raw[:2] == b"PK":
        z = zipfile.ZipFile(io.BytesIO(raw))
        xmlname = next((n for n in z.namelist() if n.lower().endswith(".xml")), None)
        if xmlname is None:
            raise InputError(f"{path.name}: no XML inside RDML archive")
        raw = z.read(xmlname)
    root = ET.fromstring(raw)

    def local(tag):
        return tag.rsplit("}", 1)[-1]

    def kids(el, name):
        return [c for c in el if local(c.tag) == name]

    def kid(el, name):
        k = kids(el, name)
        return k[0] if k else None

    def text(el, name):
        k = kid(el, name)
        return k.text.strip() if k is not None and k.text else None

    sample_type = {s.get("id"): SAMPLE_TYPES.get((text(s, "type") or "").lower(), text(s, "type") or "")
                   for s in kids(root, "sample")}
    target_dye = {}
    for t in kids(root, "target"):
        d = kid(t, "dyeId")
        target_dye[t.get("id")] = d.get("id") if d is not None else t.get("id")

    rows, cols = 8, 12
    amp = {}   # (chan, well) -> {cyc: fluor}
    melt = {}  # (chan, well) -> {tmp: fluor}
    run_names = []
    run = None
    for exp in kids(root, "experiment"):
        for r in kids(exp, "run"):
            run_names.append(r.get("id"))
            fmt = kid(r, "pcrFormat")
            if fmt is not None:
                rows = int(text(fmt, "rows") or rows)
                cols = int(text(fmt, "columns") or cols)
            if run is None:
                run = new_run(rows, cols)
            for rx in kids(r, "react"):
                rid = rx.get("id")
                if rid and rid.isdigit():
                    i = int(rid) - 1
                    well = well_name(i // cols, i % cols)
                else:
                    well = norm_well(rid or "")
                if not well:
                    continue
                s = kid(rx, "sample")
                sid = s.get("id") if s is not None else ""
                for d in kids(rx, "data"):
                    tar = kid(d, "tar")
                    tid = tar.get("id") if tar is not None else ""
                    chan = target_dye.get(tid, tid) or "Channel 1"
                    winfo = run["wells"].setdefault(well, {})
                    if sid:
                        winfo["sample"] = sid
                        if sample_type.get(sid):
                            winfo["content"] = sample_type[sid]
                    if tid:
                        winfo.setdefault("targets", set()).add(tid)
                    cq = text(d, "cq")
                    if cq:
                        try:
                            v = float(cq)
                            if v > 0:
                                run["vendor"].setdefault(chan, {}).setdefault(well, {})["cq"] = round(v, 2)
                        except ValueError:
                            pass
                    for a in kids(d, "adp"):
                        c, f = text(a, "cyc"), text(a, "fluor")
                        if c and f:
                            amp.setdefault((chan, well), {})[float(c)] = float(f)
                    for mp in kids(d, "mdp"):
                        t, f = text(mp, "tmp"), text(mp, "fluor")
                        if t and f:
                            melt.setdefault((chan, well), {})[float(t)] = float(f)
    if run is None:
        raise InputError(f"{path.name}: no runs found in RDML")
    for w in run["wells"].values():
        if "targets" in w:
            w["target"] = ", ".join(sorted(w.pop("targets")))
    _fill_curves(run, amp, melt)
    run["meta"].update(source_format="RDML export", run_name=", ".join(n for n in run_names if n))
    return run


def _fill_curves(run, amp, melt):
    """amp/melt: (chan, well) -> {x: y}. Put onto a shared x grid."""
    cyc = sorted({x for d in amp.values() for x in d})
    tmp = sorted({round(x, 2) for d in melt.values() for x in d})
    run["cycles"] = [int(x) if float(x).is_integer() else x for x in cyc]
    run["melt_temps"] = tmp
    for (chan, well), d in amp.items():
        get_channel(run, chan)["amp"][well] = [round(d[x], 2) if x in d else None for x in cyc]
    for (chan, well), d in melt.items():
        dd = {round(k, 2): v for k, v in d.items()}
        get_channel(run, chan)["melt"][well] = [round(dd[x], 2) if x in dd else None for x in tmp]


# --------------------------------------------------------------------------
# Bio-Rad data-sheet exports (CSV / XLSX)
# --------------------------------------------------------------------------
def read_table(path: Path):
    if path.suffix.lower() == ".xlsx":
        try:
            import openpyxl
        except ImportError as e:
            raise InputError("reading .xlsx exports needs openpyxl (pip install openpyxl), "
                             "or export as CSV instead") from e
        ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
        rows = [["" if v is None else str(v) for v in r] for r in ws.iter_rows(values_only=True)]
    else:
        txt = path.read_bytes().decode("utf-8-sig", errors="replace")
        rows = list(csv.reader(io.StringIO(txt)))
    # find header row (first row containing a well-like column or 'Well')
    for i, r in enumerate(rows):
        cells = [c.strip() for c in r]
        if "Well" in cells or sum(1 for c in cells if norm_well(c)) >= 2:
            return cells, [[c.strip() for c in rr] for rr in rows[i + 1:] if any(rr)]
    return [], []


def _fluor_from_name(name: str) -> str:
    m = re.search(r"_([^_]+)\.(csv|xlsx)$", name, re.I)
    return m.group(1) if m else "Channel 1"


def _num(s):
    try:
        v = float(str(s).replace(",", ""))
    except (TypeError, ValueError):
        return None
    return v if v == v and abs(v) != float("inf") else None   # drop NaN / inf


def read_biorad_exports(paths: list[Path]):
    amp, melt = {}, {}
    meta_rows = []
    peak_rows = []
    for p in paths:
        n = p.name.lower()
        if "amplification results" in n:
            kind = "amp"
        elif "melt curve rfu" in n:
            kind = "melt"
        elif "cq results" in n or "quantification summary" in n:
            kind = "meta"
        elif "melt curve peak" in n:
            kind = "peak"
        else:
            continue
        hdr, body = read_table(p)
        if not hdr:
            continue
        if kind in ("amp", "melt"):
            xcol = next((i for i, h in enumerate(hdr) if h.lower() in ("cycle", "temperature")), None)
            if xcol is None:
                continue
            wcols = [(i, norm_well(h)) for i, h in enumerate(hdr) if norm_well(h)]
            chan = _fluor_from_name(p.name)
            store = amp if kind == "amp" else melt
            for r in body:
                x = _num(r[xcol]) if xcol < len(r) else None
                if x is None:
                    continue
                for i, w in wcols:
                    y = _num(r[i]) if i < len(r) else None
                    if y is not None:
                        store.setdefault((chan, w), {})[x] = y
        else:
            idx = {h.lower(): i for i, h in enumerate(hdr)}
            for r in body:
                rec = {h: (r[i] if i < len(r) else "") for h, i in idx.items()}
                (meta_rows if kind == "meta" else peak_rows).append(rec)
    if not amp and not melt:
        raise InputError("no Bio-Rad 'Amplification Results' or 'Melt Curve RFU Results' "
                         "sheets found among the inputs")
    run = new_run()
    _fill_curves(run, amp, melt)
    for rec in meta_rows:
        w = norm_well(rec.get("well", ""))
        if not w:
            continue
        info = run["wells"].setdefault(w, {})
        for src, dst in (("sample", "sample"), ("target", "target"), ("content", "content")):
            if rec.get(src):
                info[dst] = rec[src]
        cq = _num(rec.get("cq") or rec.get("c(t)") or rec.get("ct"))
        chan = rec.get("fluor") or "Channel 1"
        if cq:
            run["vendor"].setdefault(chan, {}).setdefault(w, {})["cq"] = round(cq, 2)
    for rec in peak_rows:
        w = norm_well(rec.get("well", ""))
        tm = _num(rec.get("melt temperature"))
        if w and tm:
            chan = rec.get("fluor") or "Channel 1"
            run["vendor"].setdefault(chan, {}).setdefault(w, {}).setdefault("tm", round(tm, 2))
    run["meta"]["source_format"] = "Bio-Rad CFX data-sheet export"
    return run


# --------------------------------------------------------------------------
# Plate layout
# --------------------------------------------------------------------------
def apply_layout(run, path: Path):
    hdr, body = read_table(path)
    idx = {h.lower(): i for i, h in enumerate(hdr)}
    if "well" not in idx:
        raise InputError(f"{path.name}: layout needs a 'Well' column")
    keymap = {"sample": ["sample", "sample name", "name"], "target": ["target", "gene"],
              "content": ["content", "type"]}
    n = 0
    for r in body:
        w = norm_well(r[idx["well"]]) if idx["well"] < len(r) else None
        if not w:
            continue
        info = run["wells"].setdefault(w, {})
        for dst, opts in keymap.items():
            for o in opts:
                if o in idx and idx[o] < len(r) and r[idx[o]].strip():
                    info[dst] = r[idx[o]].strip()
                    break
        n += 1
    run["meta"]["layout_file"] = path.name
    return n


def write_layout_template(run, path: Path):
    rows, cols = run["plate"]["rows"], run["plate"]["cols"]
    with path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["Well", "Sample", "Target", "Content"])
        for r in range(rows):
            for c in range(cols):
                w = well_name(r, c)
                info = run["wells"].get(w, {})
                wr.writerow([w, info.get("sample", ""), info.get("target", ""), info.get("content", "")])


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def load(inputs: list[Path]):
    files = []
    for p in inputs:
        if p.is_dir():
            files += sorted(q for q in p.iterdir() if q.suffix.lower() in (".csv", ".xlsx", ".rdml", ".zpcr"))
        else:
            files.append(p)
    if not files:
        raise InputError("no input files found")
    pcrd = [f for f in files if f.suffix.lower() == ".pcrd"]
    if pcrd:
        raise InputError(
            f"{pcrd[0].name} is a password-protected Bio-Rad archive that only Bio-Rad software can read.\n"
            "  Open it in CFX Maestro (free download from bio-rad.com), then either\n"
            "    File > Export > Export RDML File...            -> run this script on the .rdml, or\n"
            "    Export > Export All Data Sheets (CSV format)    -> run this script on that folder.\n"
            "  Either export also carries the sample names and well types from the .pcrd.")
    zp = [f for f in files if f.suffix.lower() == ".zpcr"]
    rd = [f for f in files if f.suffix.lower() == ".rdml"]
    tb = [f for f in files if f.suffix.lower() in (".csv", ".xlsx")]
    if zp:
        run, src = read_zpcr(zp[0]), zp[0]
    elif rd:
        run, src = read_rdml(rd[0]), rd[0]
    elif tb:
        run, src = read_biorad_exports(tb), tb[0].parent if len(inputs) == 1 and inputs[0].is_dir() else tb[0]
    else:
        raise InputError("unrecognised input; expected .zpcr, .rdml, or Bio-Rad CSV/XLSX exports")
    if len(zp) + len(rd) > 1:
        print(f"note: several run files given; using {src.name}", file=sys.stderr)
    run["meta"]["source_file"] = src.name
    return run, src


def build_html(run, title: str) -> str:
    tpl = TEMPLATE.read_text(encoding="utf-8")
    payload = json.dumps(run, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    js = ANALYSIS_JS.read_text(encoding="utf-8")
    if "</script" in js.lower():
        raise InputError("analysis.js must not contain '</script' (it is inlined into the page)")
    # Fill the analysis code first, so a title or data value that happens to contain a
    # placeholder string can't pull code into the page.
    tpl = tpl.replace("__QPCR_ANALYSIS__", js)
    return tpl.replace("__QPCR_TITLE__", _html_escape(title)).replace("__QPCR_DATA__", payload)


def _html_escape(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", type=Path, help=".zpcr, .rdml, or Bio-Rad CSV/XLSX export files/folder")
    ap.add_argument("-o", "--output", type=Path, help="output .html (default: next to the input)")
    ap.add_argument("--layout", type=Path, help="CSV with Well, Sample, Target, Content columns")
    ap.add_argument("--write-layout", type=Path, metavar="CSV",
                    help="write a plate-layout CSV template for this run and exit")
    ap.add_argument("--title", help="title shown at the top of the page")
    a = ap.parse_args(argv)
    try:
        run, src = load(a.inputs)
        if a.layout:
            n = apply_layout(run, a.layout)
            print(f"layout: {n} wells from {a.layout.name}")
        if a.write_layout:
            write_layout_template(run, a.write_layout)
            print(f"wrote layout template {a.write_layout}")
            return 0
        title = a.title or src.stem
        run["meta"]["title"] = title
        run["meta"]["generated"] = dt.date.today().isoformat()
        out = a.output or (src.with_suffix(".html") if src.is_file() else src / "qpcr_viewer.html")
        out.write_text(build_html(run, title), encoding="utf-8")
    except InputError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    chans = ", ".join(c["name"] for c in run["channels"])
    print(f"read {src.name}: {len(run['cycles'])} cycles, {len(run['melt_temps'])} melt steps, channels: {chans}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
