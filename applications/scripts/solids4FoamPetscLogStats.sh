#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<'EOF'
Usage: solids4FoamPetscLogStats.sh [OPTIONS] LOG [LOG ...]

Summarise PETSc SNES/KSP metrics from one or more solids4Foam logs.

Options:
  --per-time-step  report metrics per physical time step
  -h, --help       show this help and exit

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

Per-time-step columns:
  log              input log path
  time             physical time
  fsiIterations    number of FSI coupling iteration blocks for this time step
  snesSolves       number of SNES solve reports for this time step
  totalSnes        total SNES nonlinear iterations for this time step
  maxSnes          maximum SNES nonlinear iterations in one solve
  kspSolves        number of KSP linear solve reports for this time step
  totalKsp         total KSP linear iterations for this time step
  maxKsp           maximum KSP linear iterations in one solve
  failedSnes       number of failed SNES solve reports
  retries          number of PETSc time-step retries
  resets           number of PETSc retry state resets
  divMaxIt         DIVERGED_MAX_IT reports
  divLineSearch    DIVERGED_LINE_SEARCH reports
  divFnormNan      DIVERGED_FNORM_NAN reports
  divFunctionDomain DIVERGED_FUNCTION_DOMAIN reports
  executionTime    last reported ExecutionTime for this time step, or NA
  clockTime        last reported ClockTime for this time step, or NA
EOF
}

per_time_step=false

while [[ $# -gt 0 ]]; do
    case "${1}" in
        --per-time-step)
            per_time_step=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: ${1}" >&2
            usage >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

if [[ "${per_time_step}" == true ]]; then
    :
else
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
fi

print_per_time_step_header() {
    local end_state="$1"

    printf "#end\t%s\n" "${end_state}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "#log" \
        "time" \
        "fsiIterations" \
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
}

for log_file in "$@"; do
    if [[ ! -f "${log_file}" ]]; then
        if [[ "${per_time_step}" == true ]]; then
            print_per_time_step_header "0"
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
                "${log_file}" \
                "missing" \
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                "NA" \
                "NA"
        else
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
                "${log_file}" \
                "missing" \
                0 0 0 0 0 0 0 0 0 0 0 0 0 \
                "NA" \
                "NA"
        fi
        continue
    fi

    if [[ "${per_time_step}" == true ]]; then
        if awk '/^End$/ { found = 1 } END { exit !found }' "${log_file}"; then
            print_per_time_step_header "1"
        else
            print_per_time_step_header "0"
        fi

        awk -v log_file="${log_file}" '
            function resetStep() {
                fsiIterations = 0
                snesSolves = 0
                totalSnes = 0
                maxSnes = 0
                kspSolves = 0
                totalKsp = 0
                maxKsp = 0
                failedSnes = 0
                retries = 0
                resets = 0
                divMaxIt = 0
                divLineSearch = 0
                divFnormNan = 0
                divFunctionDomain = 0
                executionTime = "NA"
                clockTime = "NA"
            }
            function flushStep() {
                if (time == "") return

                printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                    log_file,
                    time,
                    fsiIterations,
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
            /^Time =/ {
                newTime = $3
                sub(/,$/, "", newTime)

                if (time != "" && newTime != time) {
                    flushStep()
                    resetStep()
                }

                if (time == "") {
                    resetStep()
                }

                time = newTime

                if ($4 == "iteration:") {
                    fsiIterations++
                }
            }
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
            /^ExecutionTime =/ {
                executionTime = $3
                clockTime = $7
            }
            END {
                flushStep()
            }
        ' "${log_file}"
    else
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
    fi
done
