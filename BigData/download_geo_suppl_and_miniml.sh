#!/usr/bin/env bash
set -euo pipefail

# =========================
# GEO per-GSE downloader
# - Read a list of GSE IDs (one per line), download:
#   1) miniml/ (metadata)
#   2) suppl/  (supplementary data)
# - Use lftp mirror with --continue for resume
# =========================

# ---- Inputs ----
GSE_LIST_FILE="${1:-}"
DEST_ROOT="${2:-}"

if [[ -z "${GSE_LIST_FILE}" || -z "${DEST_ROOT}" ]]; then
  echo "Usage: $0 <gse_list.txt> <dest_root>"
  echo "Example: $0 gse_chunk.txt /lustre/\$USER/geo_tumor_sc"
  exit 1
fi

if [[ ! -f "${GSE_LIST_FILE}" ]]; then
  echo "ERROR: GSE list file not found: ${GSE_LIST_FILE}"
  exit 1
fi

# ---- Config ----
# GEO FTP/HTTPS root (NCBI)
REMOTE_HOST="https://ftp.ncbi.nlm.nih.gov"
REMOTE_ROOT="/geo/series"

# lftp tuning (moderate parallel to be friendly)
LFTP_PARALLEL="${LFTP_PARALLEL:-4}"
LFTP_TIMEOUT="${LFTP_TIMEOUT:-30}"
LFTP_RETRIES="${LFTP_RETRIES:-20}"

mkdir -p "${DEST_ROOT}"

echo "[INFO] GSE_LIST_FILE: ${GSE_LIST_FILE}"
echo "[INFO] DEST_ROOT:     ${DEST_ROOT}"
echo "[INFO] REMOTE_HOST:   ${REMOTE_HOST}"
echo "[INFO] LFTP_PARALLEL: ${LFTP_PARALLEL}"

# ---- Helper: compute GEO bucket (GSE20329 -> GSE20nnn) ----
gse_bucket() {
  local gse="$1"
  local num="${gse#GSE}"
  # remove any non-digits just in case
  num="$(echo "${num}" | tr -cd '0-9')"

  # 如果不是合法编号（比如表头 "GSE"），直接返回空
  [[ -z "${num}" ]] && echo "" && return 0

  local k=$(( num / 1000 ))
  # GEO: k==0 的桶叫 GSEnnn，不是 GSE0nnn
  if (( k == 0 )); then
    echo "GSEnnn"
  else
    echo "GSE${k}nnn"
  fi
}

# ---- Main loop ----
# Expect one GSE per line (e.g., GSE20329)
while read -r gse; do
  [[ -z "${gse}" ]] && continue
  # allow lines like "GSE12345\tTitle..." -> keep first field
  gse="${gse%%$'\t'*}"
  gse="$(echo "${gse}" | tr -d '[:space:]')"
  [[ -z "${gse}" ]] && continue
  [[ "${gse}" != GSE* ]] && continue
  [[ "${gse}" == "GSE" ]] && continue

  bucket="$(gse_bucket "${gse}")"
  [[ -z "${bucket}" ]] && echo "[WARN] Skip invalid line: ${gse}" && continue
  remote_dir="${REMOTE_ROOT}/${bucket}/${gse}"

  local_dest="${DEST_ROOT}/${gse}"
  mkdir -p "${local_dest}"

  echo "[INFO] Downloading ${gse}  (bucket=${bucket})"

  # We mirror two subdirs: miniml and suppl
  # -c/--continue: resume
  # --parallel: multiple simultaneous transfers within this GSE
  # NOTE: If a GSE has no suppl/ or miniml/, lftp mirror may warn but continues.

  lftp -e "
    set net:max-retries ${LFTP_RETRIES};
    set net:timeout ${LFTP_TIMEOUT};
    set xfer:clobber off;
    set ssl:verify-certificate yes;
    open ${REMOTE_HOST};

    # 1) metadata: miniml/
    mirror --continue --verbose --parallel=${LFTP_PARALLEL} \
      ${remote_dir}/miniml  ${local_dest}/miniml;

    # 2) supplementary: suppl/
    mirror --continue --verbose --parallel=${LFTP_PARALLEL} \
      ${remote_dir}/suppl   ${local_dest}/suppl;

    bye
  " || echo "[WARN] lftp failed for ${gse} (will rely on resume next run)"

done < "${GSE_LIST_FILE}"

echo "[INFO] Done. Output in: ${DEST_ROOT}"