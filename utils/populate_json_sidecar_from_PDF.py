import argparse
import shutil, subprocess
import re
import warnings
import json
from pathlib import Path
import sys

_LABEL_START = r"(?:^|(?<=\s{2}))"
_VALUE_END = r"(?=\s{2,}|\s*$)"
_FLAGS       = re.M | re.I
SI_DECIMALS = 9   # printed values carry 1-2 decimals; this keeps float noise out
PDF_FIELDS = {
    "RepetitionTime":                 dict(label="TR",               unit="ms",  scale=1e-3),
    "EchoTrainLength":                dict(label="Contrasts",        unit=None,  scale=1.0),
    "FlipAngle":                      dict(label="Flip Angle",       unit="deg", scale=1.0),
    "SliceThickness":                 dict(label="Slice Thickness",  unit="mm",  scale=1.0),
    "ParallelReductionFactorInPlane": dict(label="Accel. factor PE", unit=None,  scale=1.0),
    "PixelBandwidth":                 dict(label="Bandwidth",        unit="Hz/Px", scale=1.0),
    }
DEFAULT_FORCE = ("EchoTime",)
PRINTOUT_BLOCK = "SiemensProtocolPrintout"

class NotASiemensPrintout(ValueError):
    """Raised when the PDF does not look like a Siemens protocol printout."""

def populate_json_sidecar_from_pdf(json_path, pdf_path, output_path=None,
                                   no_strict=False, force=DEFAULT_FORCE):
    """Merge a Siemens protocol printout into a BIDS sidecar. Sidecar wins."""
    sidecar = json.loads(Path(json_path).read_text())
    text = pdf_to_text(pdf_path)
    try:
        assert_siemens(text)
    except NotASiemensPrintout:
        if not no_strict:
            raise
        warnings.warn(f"{pdf_path}: does not look like a Siemens printout; parsing anyway")

    bids, extras = parse_printout(text)

    protocol_path = extras.get("ProtocolPath")
    if protocol_path and sidecar.get("ProtocolName"):
        leaf = protocol_path.split("\\")[-1]
        if leaf != sidecar["ProtocolName"]:
            raise ValueError(
                f"PDF and sidecar describe different protocols: printout says "
                f"'{leaf}', sidecar says '{sidecar['ProtocolName']}'. "
                "Check you paired the right two files."
            )

    merged = dict(sidecar)          # copy, so the original dict is never mutated
    report = {"added": {}, "forced": {}, "confirmed": [], "conflicts": {}}

    for key, pdf_value in bids.items():
        if key not in sidecar:
            merged[key] = pdf_value
            report["added"][key] = pdf_value
        elif _same(sidecar[key], pdf_value):
            report["confirmed"].append(key)
        else:
            check_units(key, sidecar[key], pdf_value)     # raises on a x1000 error
            if key in force:
                merged[key] = pdf_value
                report["forced"][key] = {"was": sidecar[key], "now": pdf_value}
            else:
                report["conflicts"][key] = {"sidecar": sidecar[key], "pdf": pdf_value}

    merged[PRINTOUT_BLOCK] = {
        **extras,
        "SourcePDF": Path(pdf_path).name,
        "FieldsAdded": sorted(report["added"]),
        "FieldsOverwritten": sorted(report["forced"]),
    }

    if output_path:
        Path(output_path).write_text(json.dumps(merged, indent=2) + "\n")
    return merged, report

def pdf_to_text(pdf_path):
    if shutil.which("pdftotext") is None:
        raise RuntimeError("pdftotext not found on PATH (install poppler-utils)")
    out = subprocess.run(
        ["pdftotext", "-layout", str(pdf_path), "-"],
        check=True, capture_output=True, text=True,
    ).stdout
    if len(out.strip()) < 200:
        raise RuntimeError(f"{pdf_path}: no usable text layer; printout is probably image-only")
    return out

def parse_printout(text):
    """Split what the printout gives into BIDS keys and non-BIDS extras."""
    echo_times = find_echo_times(text)
    bids = parse_numeric_fields(text)
    if echo_times:
        bids["EchoTime"] = echo_times

    # 'Contrasts' and the TE count are independent statements of the same fact.
    if "EchoTrainLength" in bids and echo_times:
        if int(bids["EchoTrainLength"]) != len(echo_times):
            raise ValueError(
                f"'Contrasts' says {int(bids['EchoTrainLength'])} "
                f"but {len(echo_times)} TE values were found"
            )
        bids["EchoTrainLength"] = int(bids["EchoTrainLength"])

    extras = {
        "ProtocolPath":     find_protocol_path(text), 
        "RFPulseType":      find_text(text, "RF Pulse Type"),
        "GradientMode":     find_text(text, "Gradient Mode"),
        "PATMode":          find_text(text, "PAT mode"),
        "MultiSliceMode":   find_text(text, "Multi-slice mode"),
        "EchoSpacing":      echo_spacing(echo_times),
    }
    return bids, {k: v for k, v in extras.items() if v is not None}

def find_protocol_path(text):
    """The '\\\\USER\\...\\<protocol>' line: where this protocol lives on the scanner."""
    # [^\S\n], not \s: \s matches newlines too, so it would start the match on the
    # previous line and put m.start() on the wrong line entirely.
    m = re.search(r"^[^\S\n]*(\\\\USER\\.*)$", text, re.M)
    if not m:
        return None
    lines = text.splitlines()
    start = text[:m.start()].count("\n")
    path = lines[start].strip()
    for cont in lines[start + 1:start + 3]:
        # a wrapped continuation is short, deeply indented, and has no spaces
        if cont.strip() and " " not in cont.strip() and len(cont) - len(cont.lstrip()) > 20:
            path += cont.strip()
        else:
            break
    return path

def _same(a, b):
    """Compare sidecar and printout values, tolerating int vs float (2 vs 2.0)."""
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))
    return a == b

def find_all(text, label, unit=None):
    unit_pat = rf"\s*{re.escape(unit)}\b" if unit else ""
    pat = rf"{_LABEL_START}{label}\s{{2,}}([-\d.,]+){unit_pat}{_VALUE_END}"
    return [float(m.replace(",", ".")) for m in re.findall(pat, text, _FLAGS)]

def find_one(text, label, unit=None, scale=1.0, what=None):
    """A numeric field, which must agree everywhere the printout repeats it."""
    vals = find_all(text, label, unit)
    if not vals:
        return None
    if len(set(vals)) > 1:
        raise ValueError(
            f"{what or label}: printout disagrees with itself, found {sorted(set(vals))}"
        )
    return round(vals[0] * scale, SI_DECIMALS)

def find_text(text, label):
    """A free-text field, e.g. 'RF Pulse Type' -> 'Low SAR'. None if absent."""
    hits = re.findall(rf"{_LABEL_START}{label}\s{{2,}}(\S.*?){_VALUE_END}", text, _FLAGS)
    if not hits:
        return None
    unique = list(dict.fromkeys(h.strip() for h in hits))
    if len(unique) > 1:
        warnings.warn(f"{label}: multiple values {unique}, using the first")
    return unique[0]

def find_echo_times(text):
    pattern = rf"{_LABEL_START}TE\s+(\d+)\s{{2,}}([\d.,]+)\s*ms{_VALUE_END}"
    hits = re.findall(pattern, text, _FLAGS)
    by_index = {}
    for idx_str, val_str in hits:
        idx = int(idx_str)
        val = round(float(val_str.replace(",", ".")) * 1e-3, SI_DECIMALS)
        if idx in by_index and by_index[idx] != val:
            raise ValueError(f"TE {idx}: printout gives both {by_index[idx]} and {val} s")
        by_index[idx] = val
    if not by_index:
        return []
    missing = set(range(1, max(by_index) + 1)) - set(by_index)
    if missing:
        raise ValueError(f"echo train has gaps, missing TE {sorted(missing)}")
    return [by_index[i] for i in sorted(by_index)]

def echo_spacing(echo_times):
    """Uniform gap between echoes, or None if the train isn't uniform."""
    if len(echo_times) < 2:
        return None
    gaps = [round(b - a, SI_DECIMALS) for a, b in zip(echo_times, echo_times[1:])]
    if len(set(gaps)) > 1:
        warnings.warn(f"echo train is not uniformly spaced: gaps {gaps}")
        return None
    return gaps[0]  

def parse_numeric_fields(text):
    """Every numeric BIDS field the printout provides, already in sidecar units."""
    out = {}
    for key, spec in PDF_FIELDS.items():
        val = find_one(text, spec["label"], unit=spec["unit"],
                       scale=spec["scale"], what=key)
        if val is not None:
            out[key] = val
    return out

def check_units(key, sidecar_value, pdf_value):
    """Catch a scale mistake masquerading as a value disagreement."""
    if not (isinstance(sidecar_value, (int, float))
            and isinstance(pdf_value, (int, float))):
        return                      # lists (EchoTime) and strings aren't comparable this way
    if sidecar_value == 0:
        return                      # can't take a ratio
    ratio = pdf_value / sidecar_value
    for factor in (1e-6, 1e-3, 1e3, 1e6):
        if abs(ratio - factor) < 0.01 * factor:
            raise ValueError(
                f"{key}: sidecar has {sidecar_value}, printout gives {pdf_value} — "
                f"ratio {ratio:g} looks like a unit conversion error, not a real difference"
            )

def assert_siemens(text):
    if not re.search(r"\bSIEMENS\b", text, re.I):
        raise NotASiemensPrintout("No 'SIEMENS' banner found.")
    headings = ("Routine", "Contrast - Common", "Resolution - Common",
                "Sequence - Part 1", "Properties")
    found = sum(bool(re.search(rf"^\s*{h}\b", text, re.M)) for h in headings)
    if found < 2:
        raise NotASiemensPrintout(
            f"Siemens banner present but only {found}/{len(headings)} expected "
            "section headings found."
        )

def _brief(val):
    """Shorten long lists so the report stays one line per field."""
    if isinstance(val, list):
        return f"[{val[0]} .. {val[-1]}] ({len(val)} values)"
    return str(val)


def format_report(report):
    """A human-readable summary of what the merge did."""
    lines = []
    for key, val in report["added"].items():
        lines.append(f"  + {key:32s} {_brief(val)}")
    for key, chg in report["forced"].items():
        lines.append(f"  ! {key:32s} {_brief(chg['was'])} -> {_brief(chg['now'])}  (forced)")
    for key, chg in report["conflicts"].items():
        lines.append(f"  ? {key:32s} kept {_brief(chg['sidecar'])}, "
                     f"printout said {_brief(chg['pdf'])}")
    if report["confirmed"]:
        lines.append(f"  = confirmed: {', '.join(sorted(report['confirmed']))}")
    return "\n".join(lines) or "  (nothing to do)"
    
def parse_args():
    parser = argparse.ArgumentParser(
        description="Populate JSON sidecar from PDF"
    )
    parser.add_argument(
        "json_path",
        type=str,
        help="Path to the JSON sidecar file to populate.",
    )
    parser.add_argument(
        "pdf_path",
        type=str,
        help="Path to the PDF file to extract information from.",
    )
    parser.add_argument(
        "-o", "--output-path", type=str, default=None,
        help="Where to write the merged sidecar. Default: print to stdout.",
    )
    parser.add_argument(
        "-i", "--in-place", action="store_true",
        help="Overwrite the input sidecar in place. Cannot be combined with -o.",
    )
    parser.add_argument(
        "--force", type=str, default=",".join(DEFAULT_FORCE),
        help=f"Comma-separated keys to overwrite when the sources disagree "
             f"(default: {','.join(DEFAULT_FORCE)}; pass '' to never overwrite).",
    )
    parser.add_argument(
        "--no-strict",
        action="store_true",
        help="Warn instead of failing when the PDF is not a Siemens printout.",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if args.in_place and args.output_path:
        sys.exit("error: use either -i or -o, not both")

    destination = args.json_path if args.in_place else args.output_path
    force = tuple(k for k in args.force.split(",") if k)

    merged, report = populate_json_sidecar_from_pdf(
        args.json_path, args.pdf_path,
        output_path=destination,
        no_strict=args.no_strict,
        force=force,
    )

    if destination is None:
        print(json.dumps(merged, indent=2))
    print(format_report(report), file=sys.stderr)