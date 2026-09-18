#!/usr/bin/env python3
"""
SciFind.py
----------
Reads a BLASTP tabular output file where the 'subject sci names' column
contains N/A values, fetches the correct scientific names from NCBI Taxonomy
using the tax IDs in the 'subject tax ids' column, and writes a corrected file.

Expected BLAST fields (tab-separated, comment lines start with #):
  query acc. | subject acc. | subject tax ids | subject sci names | ...

Usage:
  python SciFind.py <input_blast_file> [output_file] --email EMAIL [--api_key KEY]

If output_file is omitted, the result is written to <input>_SciFinds.txt
"""

import argparse
import sys
import time
from pathlib import Path
from Bio import Entrez

TAX_ID_COL   = 2   # 0-based: "subject tax ids"
SCI_NAME_COL = 3   # 0-based: "subject sci names"


def fetch_sci_names(tax_ids: list[str], batch_size: int, sleep_sec: float) -> dict[str, str]:
    """Query NCBI Taxonomy for a batch of tax IDs and return {taxid: sci_name}."""
    result = {}
    for i in range(0, len(tax_ids), batch_size):
        batch = tax_ids[i : i + batch_size]
        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for rec in records:
                result[rec["TaxId"]] = rec["ScientificName"]
        except Exception as e:
            print(f"  Warning: NCBI fetch failed for batch starting at index {i}: {e}",
                  file=sys.stderr)
        time.sleep(sleep_sec)
        print(f"  Fetched {min(i + batch_size, len(tax_ids))} / {len(tax_ids)} tax IDs ...",
              end="\r", flush=True)
    print()
    return result


def collect_tax_ids(lines: list[str]) -> set[str]:
    """Extract every unique numeric tax ID from data lines."""
    ids: set[str] = set()
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= TAX_ID_COL:
            continue
        for tid in parts[TAX_ID_COL].split(";"):
            tid = tid.strip()
            if tid and tid != "N/A":
                ids.add(tid)
    return ids


def resolve_line(parts: list[str], tax_map: dict[str, str]) -> list[str]:
    """Replace the sci names field using the tax_map lookup."""
    tax_ids = [t.strip() for t in parts[TAX_ID_COL].split(";")]
    sci_names = [tax_map.get(tid, "N/A") for tid in tax_ids]
    parts[SCI_NAME_COL] = ";".join(sci_names)
    return parts


def main():
    parser = argparse.ArgumentParser(
        description="Fix missing sci names in BLAST tabular output using NCBI Taxonomy."
    )
    parser.add_argument("input", help="Input BLAST tabular file")
    parser.add_argument("output", nargs="?", help="Output file (default: <input>_SciFinds<ext>)")
    parser.add_argument("--email", required=True, help="Email address for NCBI Entrez")
    parser.add_argument("--api_key", default=None, help="NCBI API key (raises rate limit to 10 req/s)")
    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        batch_size, sleep_sec = 500, 0.11
    else:
        batch_size, sleep_sec = 200, 0.34

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_name(
        in_path.stem + "_SciFinds" + in_path.suffix
    )

    print(f"Reading {in_path} ...")
    lines = in_path.read_text().splitlines(keepends=True)

    print("Collecting tax IDs ...")
    all_tax_ids = sorted(collect_tax_ids(lines))
    print(f"  Found {len(all_tax_ids)} unique tax IDs.")

    print("Fetching scientific names from NCBI Taxonomy ...")
    tax_map = fetch_sci_names(all_tax_ids, batch_size, sleep_sec)
    print(f"  Resolved {len(tax_map)} names.")

    print("Writing corrected output ...")
    with out_path.open("w") as fh:
        for line in lines:
            if line.startswith("#") or not line.strip():
                fh.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > SCI_NAME_COL:
                parts = resolve_line(parts, tax_map)
            fh.write("\t".join(parts) + "\n")

    print(f"Done. Corrected file written to: {out_path}")


if __name__ == "__main__":
    main()
