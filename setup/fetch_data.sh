#!/usr/bin/env bash
set -euo pipefail

# Fetches the 9 files of Zenodo 5119008 into bench/data/raw/, verifies md5,
# and decompresses chrM.fa.gz into data/ref/chrM.fa for use as reference.
# Idempotent: skips files already present with correct md5.

cd "$(dirname "$0")/.."
RAW="data/raw"
REF="data/ref"
MANIFEST="data/manifest.json"
BASE="https://zenodo.org/api/records/5119008/files"

mkdir -p "$RAW" "$REF"

python3 - <<'PY' "$MANIFEST" "$RAW" "$BASE"
import hashlib, json, os, sys, urllib.request
manifest, raw, base = sys.argv[1], sys.argv[2], sys.argv[3]
with open(manifest) as f:
    m = json.load(f)
for entry in m["files"]:
    name, want_md5, want_size = entry["name"], entry["md5"], entry["size"]
    dest = os.path.join(raw, name)
    if os.path.exists(dest):
        h = hashlib.md5()
        with open(dest, "rb") as f:
            for chunk in iter(lambda: f.read(1<<16), b""):
                h.update(chunk)
        if h.hexdigest() == want_md5 and os.path.getsize(dest) == want_size:
            print(f"  ok     {name}")
            continue
        else:
            print(f"  re-dl  {name} (md5/size mismatch)")
    url = f"{base}/{name}/content"
    print(f"  fetch  {name}")
    urllib.request.urlretrieve(url, dest)
    h = hashlib.md5()
    with open(dest, "rb") as f:
        for chunk in iter(lambda: f.read(1<<16), b""):
            h.update(chunk)
    assert h.hexdigest() == want_md5, f"md5 mismatch on {name}: got {h.hexdigest()}"
    assert os.path.getsize(dest) == want_size, f"size mismatch on {name}"
print("all md5s match")
PY

if [[ ! -s "$REF/chrM.fa" ]]; then
  gunzip -c "$RAW/chrM.fa.gz" > "$REF/chrM.fa"
  echo "decompressed -> $REF/chrM.fa ($(wc -c < "$REF/chrM.fa") bytes)"
else
  echo "  ok     $REF/chrM.fa already present"
fi
