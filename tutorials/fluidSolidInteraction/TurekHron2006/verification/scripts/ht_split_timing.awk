# timestamped solids4Foam log -> per coupled step: iterations, wall, solid wall, fluid wall, KSP its/solve
{ t=$1; $1=""; line=substr($0,2) }
line ~ /^Time = [0-9.e+-]+$/ { if (step!="") out(); step=line; t0=t; nit=0; ts=0; tf=0; nk=0; sk=0 }
line ~ /^Time = .*iteration:/ { nit++ }
line ~ /^Solving the momentum equation for D/ { s0=t }
line ~ /Nonlinear solve converged|Nonlinear solve did not converge/ { ts+=t-s0 }
line ~ /^Evolving fluid model/ { f0=t }
line ~ /^Setting traction on solid interfaces/ { if (f0) tf+=t-f0; f0=0 }
line ~ /Linear solve converged/ { nk++; split(line,a," "); sk+=a[length(a)] }
function out() { printf "%-16s it=%2d wall=%6.1f solid=%6.1f fluid=%6.1f  per-it: solid %.2f fluid %.2f  ksp/solve %.0f\n", step, nit, t-t0, ts, tf, (nit?ts/(nit+1):0), (nit?tf/nit:0), (nk?sk/nk:0) }
END { out() }
