# Verify bib fields against the exact PDF each entry corresponds to.
# Mapping is manual because filenames are opaque and fuzzy matching mis-assigned.
import io, re, subprocess

BIB = r"c:/Users/afsha/Desktop/final/نهایی/مقاله latex/files/Drug_Synergy_Paper_LaTeX/references.bib"

# key -> (pdf filename or None, note)
MAP = {
    "candir2026onehot":        ("btag040.pdf",                    "published version"),
    "holbeck2017almanac":      ("3564.pdf",                       "published version"),
    "preuer2018deepsynergy":   ("bioinformatics_34_9_1538 (2).pdf","published version"),
    "kuru2022matchmaker":      ("TCBB.2021.3086702.pdf",          "published version"),
    "liu2020drugcombdb":       ("gkz1007.pdf",                    "published version"),
    "jin2025multisyn":         ("MULTI.pdf",                      "published version"),
    "ghandi2019ccle":          ("nihms-1032762.pdf",              "PMC author manuscript"),
    "theodoris2023geneformer": ("nihms-1971443.pdf",              "PMC author manuscript"),
    "janizek2018treecombo":    ("331769v1.full.pdf",              "bioRxiv (bib cites bioRxiv)"),
    "rizvi2025c2sscale":       ("C2S.pdf",                        "bioRxiv (bib cites bioRxiv)"),
    "levine2024cell2sentence": ("2023.09.11.557287.full.pdf",     "ICML camera-ready on bioRxiv"),
    "cui2024scgpt":            ("2023.04.30.538439v2.full.pdf",   "PREPRINT ONLY - bib cites Nature Methods"),
    "subramanian2017lincs":    ("136168v1.full.pdf",              "PREPRINT ONLY - bib cites Cell"),
    "wang2022deepdds":         (None, "no local PDF"),
    "elkhili2023marsy":        (None, "no local PDF - verified via Crossref"),
    "longley20035fu":          (None, "no local PDF - verified online"),
    "donawho2007abt888":       (None, "no local PDF - verified via Crossref"),
    "paul2023veliparib":       (None, "no local PDF - verified via Crossref"),
    "rodrigues20215fuorganoid":(None, "no local PDF - verified via PMC"),
    "osorio2021mito":          (None, "no local PDF - verified via publisher"),
}

def text(name, pages=2):
    r = subprocess.run(["pdftotext", "-f", "1", "-l", str(pages), "-enc", "UTF-8", name, "-"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.stdout

raw = io.open(BIB, encoding="utf-8").read()
bib = {}
for m in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", raw, re.S):
    f = {}
    for fm in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}\s*(?:,|\Z)", m.group(2), re.S):
        f[fm.group(1).lower()] = " ".join(fm.group(2).split())
    bib[m.group(1).strip()] = f

ok, warn, skip = [], [], []
for key in sorted(bib):
    f = bib[key]
    pdf, note = MAP.get(key, (None, "UNMAPPED"))
    if pdf is None:
        skip.append((key, note)); continue
    t = text(pdf)
    bad = []
    for fld in ("volume", "number", "year"):
        v = f.get(fld)
        if v and not re.search(r"(?<!\d)" + re.escape(v) + r"(?!\d)", t):
            bad.append("%s=%s" % (fld, v))
    pg = f.get("pages", "").split("--")[0].strip()
    if pg and pg.isdigit() and not re.search(r"(?<!\d)" + pg + r"(?!\d)", t):
        bad.append("firstpage=%s" % pg)
    (warn if bad else ok).append((key, pdf, note, bad))

print("=" * 96)
print("VERIFIED AGAINST THE PDF  (%d)" % len(ok))
print("=" * 96)
for k, p, n, _ in ok:
    print("  %-26s %-32s %s" % (k, p[:32], n))

print("\n" + "=" * 96)
print("FIELDS NOT FOUND ON PAGE 1-2  (%d)" % len(warn))
print("=" * 96)
for k, p, n, bad in warn:
    print("  %-26s %-32s %s" % (k, p[:32], ", ".join(bad)))
    print("  %-26s %s" % ("", n))

print("\n" + "=" * 96)
print("NO LOCAL PDF  (%d)" % len(skip))
print("=" * 96)
for k, n in skip:
    print("  %-26s %s" % (k, n))
