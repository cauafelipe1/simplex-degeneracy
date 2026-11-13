#!/usr/bin/env bash
# Downloads and decompresses NETLIB LP degeneracy benchmark instances (degen2, degen3)
# into data/mps/ using the official NETLIB emps.c utility.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
MPS_DIR="${REPO_ROOT}/data/mps"
PROCESSED_DIR="${REPO_ROOT}/data/processed"
NETLIB_BASE_URL="https://www.netlib.org/lp/data"

mkdir -p "${MPS_DIR}" "${PROCESSED_DIR}"

for cmd in curl gcc; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Error: '${cmd}' is required to download and decompress NETLIB instances." >&2
        exit 1
    fi
done

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "[1/3] Downloading official NETLIB decompressor (emps.c)..."
curl -fsSL "${NETLIB_BASE_URL}/emps.c" -o "${TMP_DIR}/emps.c"
gcc -O2 "${TMP_DIR}/emps.c" -o "${TMP_DIR}/emps"

for inst in degen2 degen3; do
    echo "[2/3] Downloading and decompressing ${inst}.mps..."
    curl -fsSL "${NETLIB_BASE_URL}/${inst}" -o "${TMP_DIR}/${inst}.emps"
    "${TMP_DIR}/emps" "${TMP_DIR}/${inst}.emps" > "${MPS_DIR}/${inst}.mps"
done

echo "[3/3] Converting MPS files to sparse Octave .mat format..."
if command -v octave >/dev/null 2>&1; then
    octave --no-gui --eval "
        addpath('${REPO_ROOT}/src');
        mps_parser('${MPS_DIR}/degen2.mps', '${PROCESSED_DIR}/degen2.mat', true);
        mps_parser('${MPS_DIR}/degen3.mps', '${PROCESSED_DIR}/degen3.mat', true);
    "
else
    echo "Warning: 'octave' not found in PATH. Raw .mps files saved in ${MPS_DIR};"
    echo "         they will be converted automatically when running the benchmark."
fi

echo "Done! NETLIB instances are ready in ${MPS_DIR}."
