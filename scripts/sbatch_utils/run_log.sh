#!/usr/bin/env bash
# Central run-log helpers, sourced from sbatch scripts.
#
# Required before start_run_log:
#   KEY, INPUT_PATH, RESULT_DIR, CMD
# Optional:
#   METHOD (default: basename of RESULT_DIR)
#   STAINING, CP_VERSION, CONFIDENCE, PARAMS (default NA)
# Override RUN_LOG via env if needed.

RUN_LOG="${RUN_LOG:-/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/run_log.tsv}"
LOCK_FILE="${RUN_LOG}.lock"
RUN_LOG_HEADER=$'start_iso\tend_iso\telapsed_s\trc\tjobid\tjobname\tkey\tmethod\tcp_version\tstaining\tconfidence\tparams\tinput_path\tresult_dir\thost\tnodelist\tsubmit_dir\tcmd'

start_run_log() {
  : "${KEY:?}" "${INPUT_PATH:?}" "${RESULT_DIR:?}" "${CMD:?}"
  METHOD="${METHOD:-$(basename "${RESULT_DIR}")}"
  STAINING="${STAINING:-NA}"
  CP_VERSION="${CP_VERSION:-NA}"
  CONFIDENCE="${CONFIDENCE:-NA}"
  PARAMS="${PARAMS:-NA}"
  JOBID="${SLURM_JOB_ID:-NA}"
  JOBNAME="${SLURM_JOB_NAME:-NA}"
  NODELIST="${SLURM_JOB_NODELIST:-NA}"
  SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
  HOST="$(hostname -f 2>/dev/null || hostname)"
  START_ISO="$(date -Is)"
  START_EPOCH="$(date +%s)"
  mkdir -p "$(dirname "${RUN_LOG}")"
  trap '_write_run_log $?' EXIT
}

_write_run_log() {
  local rc="$1" end_iso elapsed_s hdr
  end_iso="$(date -Is)"
  elapsed_s=$(( $(date +%s) - START_EPOCH ))
  (
    flock -x 200
    if [ ! -s "${RUN_LOG}" ]; then
      printf '%s\n' "${RUN_LOG_HEADER}" >> "${RUN_LOG}"
    else
      IFS= read -r hdr < "${RUN_LOG}" || true
      [ "${hdr}" = "${RUN_LOG_HEADER}" ] || \
        printf 'run_log: header mismatch in %s\n' "${RUN_LOG}" >&2
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${START_ISO}" "${end_iso}" "${elapsed_s}" "${rc}" \
      "${JOBID}" "${JOBNAME}" "${KEY}" "${METHOD}" "${CP_VERSION}" "${STAINING}" "${CONFIDENCE}" "${PARAMS}" \
      "${INPUT_PATH}" "${RESULT_DIR}" "${HOST}" "${NODELIST}" "${SUBMIT_DIR}" "${CMD}" >> "${RUN_LOG}"
  ) 200>>"${LOCK_FILE}"
}