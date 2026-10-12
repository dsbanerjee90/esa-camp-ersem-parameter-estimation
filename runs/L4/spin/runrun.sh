#!/usr/bin/env bash
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
set -euo pipefail

GOTM_BIN="${GOTM_BIN:-eat-gotm}"

echo "================================================================================"
echo "CAMP_ARGO32 SPIN-UP"
echo "Working directory: $(pwd)"
echo "GOTM executable:   ${GOTM_BIN}"
echo "================================================================================"

for required_file in gotm.yaml fabm.yaml; do
    if [[ ! -s "${required_file}" ]]; then
        echo "ERROR: missing or empty file: ${required_file}" >&2
        exit 1
    fi
done

if ! command -v "${GOTM_BIN}" >/dev/null 2>&1; then
    echo "ERROR: executable not found in PATH: ${GOTM_BIN}" >&2
    exit 1
fi

python - <<'PY'
from pathlib import Path
import re
import yaml

path = Path("gotm.yaml")
text = path.read_text(encoding="utf-8")

# GOTM uses the literal string "off". Prevent PyYAML from treating
# off/on/yes/no as YAML 1.1 Boolean values during validation.
class GotmLoader(yaml.SafeLoader):
    pass

GotmLoader.yaml_implicit_resolvers = {
    key: list(value)
    for key, value in yaml.SafeLoader.yaml_implicit_resolvers.items()
}

boolean_tag = "tag:yaml.org,2002:bool"

for first_character, entries in list(
    GotmLoader.yaml_implicit_resolvers.items()
):
    GotmLoader.yaml_implicit_resolvers[first_character] = [
        (tag, pattern)
        for tag, pattern in entries
        if tag != boolean_tag
    ]

GotmLoader.add_implicit_resolver(
    boolean_tag,
    re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"),
    list("tTfF"),
)

cfg = yaml.load(
    text,
    Loader=GotmLoader,
) or {}

time_cfg = cfg.get("time", {})
restart_cfg = cfg.get("restart", {})
fabm_cfg = cfg.get("fabm", {})

missing = [
    key
    for key in ("start", "stop")
    if not time_cfg.get(key)
]

if missing:
    raise SystemExit(
        "ERROR: missing required GOTM time setting(s): "
        + ", ".join(missing)
    )

if restart_cfg != {"load": False}:
    raise SystemExit(
        "ERROR: spin gotm.yaml restart section must contain "
        "only: load: false"
    )

if not fabm_cfg.get("use"):
    raise SystemExit(
        "ERROR: FABM is not enabled in gotm.yaml"
    )

print("Configuration checks passed")
print("  Start:", time_cfg["start"])
print("  Stop: ", time_cfg["stop"])
print("  FABM:  fabm.yaml")
print("  Restart input: disabled")
PY

if [[ -e restart.nc ]]; then
    echo "ERROR: restart.nc already exists." >&2
    echo "Refusing to overwrite a previous spin result." >&2
    exit 1
fi

echo
echo "Starting free GOTM spin-up..."
"${GOTM_BIN}"

if [[ ! -s restart.nc ]]; then
    echo "ERROR: spin completed without producing restart.nc" >&2
    exit 1
fi

echo
echo "SPIN-UP COMPLETED"
echo "Restart created: $(pwd)/restart.nc"
ls -lh restart.nc
