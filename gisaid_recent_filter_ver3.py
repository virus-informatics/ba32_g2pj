#!/usr/bin/env python
# coding: utf-8

import lzma
import sys
import logging
from datetime import date

argvs = sys.argv

# Usage: python3 gisaid_recent_filter_ver4.py <path_to_sequences.tar.xz> [log_file]
full_data = lzma.open(argvs[1])
log_path = argvs[2] if len(argvs) > 2 else "gisaid_filter.log"

# ---- Fixed collection-date window ----
start_date = date.fromisoformat("2025-08-01")
end_date   = date.fromisoformat("2026-07-01")

# ---- Header format ----
# >hCoV-19/Country/ID/Year|COLLECTION_DATE|SUBMISSION_DATE
COLLECTION_DATE_FIELD = 1

# ---- Logging setup: warnings go to a file, NOT the terminal ----
logging.basicConfig(
    filename=log_path,
    filemode="w",
    level=logging.INFO,
    format="%(message)s",
)
logger = logging.getLogger(__name__)

PROGRESS_EVERY = 200_000  # print a heartbeat to stderr every N headers, so you know it's alive

n_total = 0
n_kept = 0
n_bad_date = 0
n_bad_date_relevant = 0  # partial dates whose year overlaps the target window (worth a closer look)

target_years = set(range(start_date.year, end_date.year + 1))

keep = 'no'
for line in full_data:
    try:
        line = line.decode()
    except ValueError:
        continue

    if line.startswith('>'):
        n_total += 1
        if n_total % PROGRESS_EVERY == 0:
            print(f"...processed {n_total:,} headers so far (kept {n_kept:,})", file=sys.stderr)

        fields = line.strip().split('|')
        seq_date = None
        raw_field = fields[COLLECTION_DATE_FIELD] if len(fields) > COLLECTION_DATE_FIELD else ""

        if len(fields) > COLLECTION_DATE_FIELD:
            try:
                seq_date = date.fromisoformat(raw_field)
            except ValueError:
                seq_date = None

        if seq_date is None:
            n_bad_date += 1
            # Only log (and only count as "relevant") if the partial date's year
            # falls within the target window's year range -- e.g. "2025-08" or "2026"
            # -- since that's the only case where missing precision could have mattered.
            year_str = raw_field[:4]
            if year_str.isdigit() and int(year_str) in target_years:
                n_bad_date_relevant += 1
                logger.info(f"RELEVANT partial/bad date, excluded: {line.strip()}")
            keep = 'no'
        elif start_date <= seq_date <= end_date:
            keep = 'yes'
            n_kept += 1
            print(fields[0], end="\n")
        else:
            keep = 'no'
    else:
        if keep == 'yes':
            print(line, end="")

summary = (
    f"Total headers seen: {n_total:,}\n"
    f"Kept ({start_date} to {end_date}): {n_kept:,}\n"
    f"Headers with unparseable/partial date: {n_bad_date:,}\n"
    f"  of which year overlapped target window (worth reviewing): {n_bad_date_relevant:,}\n"
    f"Full list of relevant excluded headers written to: {log_path}\n"
)
print(summary, file=sys.stderr)
logger.info("\n" + summary)
