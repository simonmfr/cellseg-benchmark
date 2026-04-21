#!/usr/bin/env bash
# Central run-log helpers, sourced from sbatch scripts.
#
# Required before start_run_log:
#   KEY, METHOD, INPUT_PATH, RESULT_DIR, CMD
# Optional (default NA):
#   STAINING, CP_VERSION, CONFIDENCE, PARAMS
# Override RUN_LOG via env if needed.

RUN_LOG="${RUN_LOG:-/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/job_runs.tsv}"
LOCK_FILE="${RUN_LOG}.lock"

start_run_log() {
  : "${KEY:?}" "${METHOD:?}" "${INPUT_PATH:?}" "${RESULT_DIR:?}" "${CMD:?}"
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
  local rc="$1" end_iso elapsed_s
  end_iso="$(date -Is)"
  elapsed_s=$(( $(date +%s) - START_EPOCH ))
  (
    flock -x 200
    if [ ! -s "${RUN_LOG}" ]; then
      printf "start_iso\tend_iso\telapsed_s\trc\tjobid\tjobname\tkey\tmethod\tcp_version\tstaining\tconfidence\tparams\tinput_path\tresult_dir\thost\tnodelist\tsubmit_dir\tcmd\n" >> "${RUN_LOG}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${START_ISO}" "${end_iso}" "${elapsed_s}" "${rc}" \
      "${JOBID}" "${JOBNAME}" "${KEY}" "${METHOD}" "${CP_VERSION}" "${STAINING}" "${CONFIDENCE}" "${PARAMS}" \
      "${INPUT_PATH}" "${RESULT_DIR}" "${HOST}" "${NODELIST}" "${SUBMIT_DIR}" "${CMD}"
  ) 200>>"${LOCK_FILE}" >> "${RUN_LOG}"
}