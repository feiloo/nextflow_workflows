#!/usr/bin/env python3

# this script was written by an AI

#!/usr/bin/env python3
import sys
import re

def parse_bibtex_entries(text):
    pattern = r'@(\w+)\s*\{\s*([^,]+),(.+?)\n\}'
    matches = re.findall(pattern, text, flags=re.S)
    entries = []

    for entry_type, key, body in matches:
        fields = {}
        parts = re.split(r',(?![^{]*\})', body)
        for part in parts:
            if '=' in part:
                k, v = part.split('=', 1)
                k = k.strip().lower()
                v = v.strip().rstrip(',').strip()
                v = re.sub(r'^[{\"]|[}\"]$', '', v).strip()
                fields[k] = v
        fields["entry_type"] = entry_type.lower()
        fields["cite_key"] = key.strip()
        entries.append(fields)

    return entries


# words commonly belonging to multi-word last names
LASTNAME_PREFIXES = {
    "van", "von", "der", "de", "da", "del", "della",
    "di", "la", "le", "du", "dos", "das", "el"
}


def format_author(name):
    """
    Return an author in 'First Middle Last' form (full names, no initials).
    Handles:
      - First Middle Last
      - Last, First Middle
      - multi-word last names ('van der Waals')
    """
    name = name.strip()

    # Case 1: "Last, First Middle"
    if ',' in name:
        last, firsts = [p.strip() for p in name.split(',', 1)]
        return f"{firsts} {last}".strip()

    # Case 2: "First Middle Last" (possibly multi-word last name)
    parts = name.split()
    if len(parts) == 1:
        return parts[0]

    # Detect last name (preserve multiword last names)
    last_start = len(parts) - 1
    while last_start > 0 and parts[last_start - 1].lower() in LASTNAME_PREFIXES:
        last_start -= 1

    firsts = " ".join(parts[:last_start])
    last = " ".join(parts[last_start:])
    return f"{firsts} {last}".strip()


def format_ieee(entry):
    authors_raw = entry.get("author", "").replace("\n", " ")
    authors = [a.strip() for a in authors_raw.split(" and ") if a.strip()]
    authors_fmt = ", ".join(format_author(a) for a in authors)

    title = entry.get("title", "")
    journal = entry.get("journal") or entry.get("booktitle") or ""
    year = entry.get("year", "")
    volume = entry.get("volume")
    number = entry.get("number")
    pages = entry.get("pages")
    publisher = entry.get("publisher", "") if entry["entry_type"] == "book" else ""

    parts = [authors_fmt, f'"{title}"']
    if journal: parts.append(journal)
    if publisher: parts.append(publisher)
    if volume: parts.append(f"vol. {volume}")
    if number: parts.append(f"no. {number}")
    if pages: parts.append(f"pp. {pages}")
    if year: parts.append(year)

    return ", ".join(p for p in parts if p)


def main():
    if len(sys.argv) != 2:
        print("Usage: python bib2ieee.py file.bib")
        sys.exit(1)

    with open(sys.argv[1], "r", encoding="utf-8") as f:
        text = f.read()

    entries = parse_bibtex_entries(text)
    for e in entries:
        print(format_ieee(e))


if __name__ == "__main__":
    main()

