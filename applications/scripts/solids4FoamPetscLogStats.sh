#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<'EOF'
Usage: solids4FoamPetscLogStats.sh LOG [LOG ...]

Summarise PETSc SNES/KSP metrics from one or more solids4Foam logs.

Columns:
  log              input log path
  end              yes if the solver log reached End
  snesSolves       number of SNES solve reports
  totalSnes        total SNES nonlinear iterations
  maxSnes          maximum SNES nonlinear iterations in one solve
  kspSolves        number of KSP linear solve reports
  totalKsp         total KSP linear iterations
  maxKsp           maximum KSP linear iterations in one solve
  failedSnes       number of failed SNES solve reports
  retries          number of PETSc time-step retries
  resets           number of PETSc retry state resets
  divMaxIt         DIVERGED_MAX_IT reports
  divLineSearch    DIVERGED_LINE_SEARCH reports
  divFnormNan      DIVERGED_FNORM_NAN reports
  divFunctionDomain DIVERGED_FUNCTION_DOMAIN reports
  executionTime    final reported ExecutionTime, or NA
  clockTime        final reported ClockTime, or NA
EOF
}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "log" \
    "end" \
    "snesSolves" \
    "totalSnes" \
    "maxSnes" \
    "kspSolves" \
    "totalKsp" \
    "maxKsp" \
    "failedSnes" \
    "retries" \
    "resets" \
    "divMaxIt" \
    "divLineSearch" \
    "divFnormNan" \
    "divFunctionDomain" \
    "executionTime" \
    "clockTime"

for log_file in "$@"; do
    if [[ ! -f "${log_file}" ]]; then
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "${log_file}" \
            "missing" \
            0 0 0 0 0 0 0 0 0 0 0 0 0 \
            "NA" \
            "NA"
        continue
    fi

    awk -v log_file="${log_file}" '
        /Linear solve .* iterations/ {
            totalKsp += $NF
            if ($NF > maxKsp) maxKsp = $NF
            kspSolves++
        }
        /Nonlinear solve (converged|did not converge).* iterations/ {
            totalSnes += $NF
            if ($NF > maxSnes) maxSnes = $NF
            snesSolves++
        }
        /Nonlinear solve did not converge/ {
            failedSnes++
        }
        /Retrying the failed PETSc time step/ {
            retries++
        }
        /Resetting PETSc SNES\/KSP retry state/ {
            resets++
        }
        /DIVERGED_MAX_IT/ {
            divMaxIt++
        }
        /DIVERGED_LINE_SEARCH/ {
            divLineSearch++
        }
        /DIVERGED_FNORM_NAN/ {
            divFnormNan++
        }
        /DIVERGED_FUNCTION_DOMAIN/ {
            divFunctionDomain++
        }
        /^End$/ {
            reachedEnd = 1
        }
        /^ExecutionTime =/ {
            executionTime = $3
            clockTime = $7
        }
        END {
            endState = "no"
            if (reachedEnd) endState = "yes"
            if (executionTime == "") executionTime = "NA"
            if (clockTime == "") clockTime = "NA"

            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                log_file,
                endState,
                snesSolves,
                totalSnes,
                maxSnes,
                kspSolves,
                totalKsp,
                maxKsp,
                failedSnes,
                retries,
                resets,
                divMaxIt,
                divLineSearch,
                divFnormNan,
                divFunctionDomain,
                executionTime,
                clockTime
        }
    ' "${log_file}"
done
