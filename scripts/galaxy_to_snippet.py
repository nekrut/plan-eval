#!/usr/bin/env python3
"""
Deterministic Cheetah-XML → bash snippet extractor for Galaxy IUC tool wrappers.

Pipeline:
 1. Parse the main XML, find <macros><import>FILE</import></macros>.
 2. Build a token table from those macro files (<token name="@NAME@">value</token>).
 3. Substitute @TOKEN@ references in the <command> block, recursively (up to N
    expansions to avoid cycles).
 4. Strip Cheetah conditionals (#if/#elif/#else/#end if), #set/#silent/#for, and
    ## comment lines. Bias to minimal-invocation: never expand conditional bodies.
 5. Drop pure shell-glue lines (ln -s, mkdir, cd, &&, ||, error handlers).
 6. Apply project-specific variable substitutions (Galaxy filenames → our paths).
 7. Drop residual unbound $var references.
 8. Collapse whitespace.

Usage:
  galaxy_to_snippet.py <main.xml> [<macros_dir>...]

If --macros-dir args are provided, the script searches them for the imported
macro filenames (in addition to the main XML's directory).
"""
from __future__ import annotations
import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Project-specific variable substitution (longest patterns first).
SUBS: list[tuple[str, str]] = [
    (r"\\?\$\{GALAXY_SLOTS:-\d+\}", "4"),
    (r"\\?\$\{GALAXY_MEMORY_MB[^}]*\}", "4096"),
    (r"'\$reference_fasta_fn'", "data/ref/chrM.fa"),
    (r"'\$ref_file'",            "data/ref/chrM.fa"),
    (r"'\$reference'",           "data/ref/chrM.fa"),
    (r"'\$input1'",              "results/{sample}.bam"),
    (r"'\$input'",               "results/{sample}.bam"),
    (r"'\$reads'",               "results/{sample}.bam"),
    (r"'\$bam_input'",           "results/{sample}.bam"),
    (r"\breads\.bam\b",          "results/{sample}.bam"),
    (r"\bvariants\.vcf\b",       "results/{sample}.vcf"),
    (r"'\$output1'",             "results/{sample}.out"),
    (r"'\$output'",              "results/{sample}.out"),
    # Strip Galaxy macro references that didn't expand
    (r"@[A-Z_][A-Z0-9_]*@",      ""),
    # Strip backslash-escaped $var (after named subs handled the rest)
    (r"\\\$([A-Za-z_][A-Za-z0-9_]*)", r"\1"),
]

CHEETAH_VAR_RE = re.compile(r"(?<!\\)\$[A-Za-z_][A-Za-z0-9_.]*")

DROP_PREFIXES = (
    "ln -s", "ln -sf", "mkdir ", "cd ", "rm ", "cp ",
    "&&", "||", "2>&1",
    "tool_exit_code", "exit ", "cat ", "echo ",
    "## ", "#set ", "#silent ", "#for ", "#end for", "#end if",
    "#if ", "#elif ", "#else", "#try", "#end try",
)
GLUE_TOKEN_RE = re.compile(r"\s*&&\s*$|\s*\|\|\s*$|\s*2>&1\s*$")


def _safe_parse(xml_path: Path) -> ET.Element | None:
    """Parse an XML file, tolerating <macros>'s unbalanced root if it's a
    macro file rather than a tool wrapper."""
    try:
        return ET.parse(xml_path).getroot()
    except ET.ParseError:
        # Some macro files use <macros> as the implicit root; try wrapping
        text = xml_path.read_text()
        try:
            return ET.fromstring(f"<wrapper>{text}</wrapper>")
        except ET.ParseError:
            return None


def load_tokens(xml_paths: list[Path]) -> dict[str, str]:
    """Build {token_name: value} from one or more macro XML files."""
    tokens: dict[str, str] = {}
    for path in xml_paths:
        root = _safe_parse(path)
        if root is None:
            continue
        for tok in root.iter("token"):
            name = tok.get("name", "")
            value = (tok.text or "").strip()
            if name:
                tokens[name] = value
    return tokens


def find_macro_imports(xml_path: Path, search_dirs: list[Path]) -> list[Path]:
    """Resolve <macros><import>...</import></macros> filenames in search_dirs."""
    root = _safe_parse(xml_path)
    if root is None:
        return []
    found: list[Path] = []
    for imp in root.iter("import"):
        fname = (imp.text or "").strip()
        if not fname:
            continue
        for d in search_dirs:
            candidate = d / fname
            if candidate.exists():
                found.append(candidate)
                break
    return found


def expand_tokens(s: str, tokens: dict[str, str], max_passes: int = 5) -> str:
    """Repeatedly substitute @NAME@ → token value until fixed point or max_passes."""
    for _ in range(max_passes):
        replaced = False
        for name, value in tokens.items():
            if name in s:
                s = s.replace(name, value)
                replaced = True
        if not replaced:
            break
    return s


def extract_command(xml_path: Path) -> str:
    root = _safe_parse(xml_path)
    if root is None:
        raise RuntimeError(f"{xml_path}: failed to parse")
    cmd = root.find(".//command")
    if cmd is None or cmd.text is None:
        raise RuntimeError(f"{xml_path}: no <command> element with text")
    return cmd.text


def strip_cheetah_conditionals(s: str) -> str:
    lines = s.splitlines()
    out, depth = [], 0
    for line in lines:
        st = line.strip()
        if st.startswith("#if "):
            depth += 1
            continue
        if st == "#end if":
            depth = max(0, depth - 1)
            continue
        if depth == 0:
            out.append(line)
    return "\n".join(out)


def drop_glue_lines(s: str) -> str:
    out = []
    for line in s.splitlines():
        st = line.strip()
        if not st:
            continue
        if any(st.startswith(p) for p in DROP_PREFIXES):
            continue
        if st.startswith("|| ("):
            continue
        line = GLUE_TOKEN_RE.sub("", line)
        if line.strip():
            out.append(line.rstrip())
    return "\n".join(out)


def apply_subs(s: str) -> str:
    for pat, rep in SUBS:
        s = re.sub(pat, rep, s)
    s = CHEETAH_VAR_RE.sub("", s)
    return s


def collapse(s: str) -> str:
    s = re.sub(r"\\\n\s*", " ", s)
    out_lines = [re.sub(r"\s+", " ", l).strip() for l in s.splitlines()]
    out_lines = [l for l in out_lines if l]
    return "\n".join(out_lines)


# Lines we drop in a final pass after collapse: residual Cheetah `#def`, leading
# `@MACRO@` references, lone parens, and lines that start with `${` (unbound
# Cheetah var reference).
_NOISE_LINE_RE = re.compile(r"^(#|@|\(|\$\{|##|//)")


def post_clean(s: str) -> str:
    out = []
    for line in s.splitlines():
        if _NOISE_LINE_RE.match(line.strip()):
            continue
        # Strip stray empty `${...}` and `$.foo` residue
        cleaned = re.sub(r"\$\{[^}]*\}", "", line)
        cleaned = re.sub(r"\s{2,}", " ", cleaned).strip()
        if cleaned:
            out.append(cleaned)
    return "\n".join(out)


def extract_snippet(xml_path: Path, search_dirs: list[Path]) -> str:
    macro_paths = find_macro_imports(xml_path, search_dirs)
    tokens = load_tokens(macro_paths)
    raw = extract_command(xml_path)
    s = expand_tokens(raw, tokens)
    s = strip_cheetah_conditionals(s)
    s = drop_glue_lines(s)
    s = apply_subs(s)
    s = collapse(s)
    s = post_clean(s)
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("xml")
    ap.add_argument("--macros-dir", action="append", default=[],
                    help="Additional directory to search for macro imports (repeatable)")
    args = ap.parse_args()
    xml_path = Path(args.xml)
    search_dirs = [xml_path.parent] + [Path(d) for d in args.macros_dir]
    print(extract_snippet(xml_path, search_dirs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
