"""Shared helpers for the functionalhomology pipeline scripts."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable

_ACC_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")


def open_text(path: str | Path):
    """Open a plain or gzip-compressed text file."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def accession_from_header(header: str) -> str:
    """Extract a UniProt accession from a FASTA header.

    Handles ``sp|P12345|NAME_ECOLI ...`` and ``tr|A0A123|...``.  Falls back to
    the first whitespace-delimited token, stripping a ``>`` if present.
    """
    header = header.lstrip(">").strip()
    m = _ACC_RE.match(header)
    if m:
        return m.group(1)
    return header.split()[0] if header else ""


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file as (accession, sequence) pairs.

    Raises a clear error on empty files or duplicate accessions.
    """
    records: list[tuple[str, str]] = []
    seen: set[str] = set()
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal header, chunks
        if header is None:
            return
        acc = accession_from_header(header)
        if not acc:
            raise ValueError(f"Could not parse accession from header: {header!r}")
        if acc in seen:
            raise ValueError(f"Duplicate accession {acc!r} in {path}")
        seen.add(acc)
        records.append((acc, "".join(chunks).replace(" ", "").replace("\t", "")))
        header, chunks = None, []

    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line
            else:
                chunks.append(line.strip())
    flush()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def fasta_names(path: str | Path) -> dict[str, str]:
    """Map accession -> human-readable name/description from FASTA headers."""
    names: dict[str, str] = {}
    with open_text(path) as fh:
        header: str | None = None
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    _store_name(header, names)
                header = line
        if header is not None:
            _store_name(header, names)
    return names


def _store_name(header: str, names: dict[str, str]) -> None:
    acc = accession_from_header(header)
    parts = header.lstrip(">").split(maxsplit=1)
    entry = parts[0].split("|")[-1] if "|" in parts[0] else parts[0]
    desc = parts[1] if len(parts) > 1 else ""
    names[acc] = f"{entry} — {desc}" if desc else entry


def parse_obo(path: str | Path):
    """Parse an OBO file.

    Returns ``(blocks, parents)`` where ``blocks`` maps a GO id to its raw
    ``[Term]`` lines and ``parents`` maps a GO id to the set of GO ids linked
    by ``is_a`` or ``relationship: part_of``.
    """
    blocks: dict[str, list[str]] = {}
    parents: dict[str, set[str]] = {}
    current: list[str] | None = None
    current_id: str | None = None
    in_term = False

    def flush() -> None:
        nonlocal current, current_id
        if in_term and current_id:
            blocks[current_id] = current or []
        current, current_id = None, None

    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("[Term]"):
                flush()
                in_term = True
                current = [line]
                continue
            if line.startswith("[") and line != "[Term]":
                flush()
                in_term = False
                continue
            if not in_term:
                continue
            current = current or []
            current.append(line)
            if line.startswith("id: "):
                current_id = line[4:].strip()
                parents.setdefault(current_id, set())
            elif line.startswith("is_a: "):
                parent = line[6:].split("!", 1)[0].strip()
                if current_id and parent.startswith("GO:"):
                    parents.setdefault(current_id, set()).add(parent)
            elif line.startswith("relationship: part_of "):
                parent = line.split("part_of", 1)[1].split("!", 1)[0].strip()
                if current_id and parent.startswith("GO:"):
                    parents.setdefault(current_id, set()).add(parent)
    flush()
    return blocks, parents


def obo_names(blocks: dict[str, list[str]]) -> dict[str, str]:
    """Extract GO id -> name from parsed OBO blocks."""
    names: dict[str, str] = {}
    for go_id, block in blocks.items():
        for line in block:
            if line.startswith("name: "):
                names[go_id] = line[6:].strip()
                break
    return names


def ancestor_closure(terms: Iterable[str], parents: dict[str, set[str]]) -> set[str]:
    """Return terms plus all ancestors reachable via is_a/part_of."""
    seen: set[str] = set()
    stack = list(terms)
    while stack:
        term = stack.pop()
        if term in seen:
            continue
        seen.add(term)
        stack.extend(parents.get(term, ()))
    return seen
