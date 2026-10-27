# Final-cycle comparison of the coupled closed-loop pipe runs with the
# monolithic R-L-surrogate reference and with the original runs.
#   julia plotPipeClosedLoop.jl [runs] [outdir]
# Lines: this code. Hollow circles: the original run (every 20th step).
# Black dashed: the monolithic reference (pipe replaced by its R-L surrogate).
ENV["GKSwstype"] = "100"
using Plots

here = @__DIR__
runs = length(ARGS) >= 1 ? ARGS[1] : joinpath(here, "runs")
out = length(ARGS) >= 2 ? ARGS[2] : joinpath(here, "plots")
mkpath(out)
refdir = joinpath(here, "..", "..", "verification", "reference")

function readtable(path, key)
    names = String[]; rows = Vector{Vector{Float64}}()
    for l in eachline(path)
        if startswith(l, "#")
            p = split(l[2:end]); if !isempty(p) && p[1] == key; names = String.(p); end
        elseif !isempty(strip(l)); push!(rows, parse.(Float64, split(l))); end
    end
    return names, reduce(hcat, rows)'
end

mono = let ls = readlines(joinpath(refdir, "closedLoop_monolithicRK4_dt6p25e-6.csv"))
    n = String.(split(ls[1], ",")); M = reduce(hcat, [parse.(Float64, split(l, ",")) for l in ls[2:end]])'
    (n, M)
end
mcol(n) = mono[2][:, findfirst(==(n), mono[1])]
tm = mcol("time"); fin_m = tm .>= 1.6

function load(scheme, order, dt)
    f = joinpath(runs, scheme, string(order), dt, "postProcessing", "circulationCoupling", "0", "coupling.dat")
    isfile(f) || return nothing
    readtable(f, "time")
end
function loadold(scheme, order, dt)
    nm = scheme * ["First", "Second", "Third", "Fourth"][order] * "Order"
    f = joinpath(refdir, "pipeClosedLoop", nm, "dt_" * dt, "stateLog.dat")
    isfile(f) || return nothing
    readtable(f, "time")
end

cases = [("implicit", 2, "5e-3", :darkorange), ("explicit", 4, "5e-3", :crimson), ("explicit", 2, "5e-3", :royalblue)]

kw = (framestyle=:box, legend=:topright, xlabel="t (s)")
pPo = plot(; ylabel="Po (Pa)", title="(a) pipe outlet pressure", kw...)
pPi = plot(; ylabel="Pi (Pa)", title="(b) pipe inlet pressure", kw...)
pQ  = plot(; ylabel="Q (m3/s)", title="(c) flow through the pipe", kw...)
pPV = plot(; xlabel="Vv (mL)", ylabel="Pv (kPa)", title="(d) ventricular PV loop", framestyle=:box, legend=:topright)

plot!(pPo, tm[fin_m], mcol("Po")[fin_m]; color=:black, ls=:dash, lw=2, label="monolithic R-L reference")
plot!(pPi, tm[fin_m], mcol("Pi")[fin_m]; color=:black, ls=:dash, lw=2, label="")
plot!(pQ, tm[fin_m], mcol("Qpipe")[fin_m]; color=:black, ls=:dash, lw=2, label="")

for (sch, o, dt, c) in cases
    d = load(sch, o, dt); d === nothing && continue
    n, M = d; col(x) = M[:, findfirst(==(x), n)]; t = col("time"); f = t .>= 1.6
    lab = "$sch $o, dt=$dt"
    plot!(pPo, t[f], col("C1.P")[f]; color=c, lw=2, label=lab)
    plot!(pPi, t[f], col("C4.P")[f]; color=c, lw=2, label="")
    plot!(pQ, t[f], col("inlet.Q")[f]; color=c, lw=2, label="")
    plot!(pPV, col("ventricle.V")[f]*1e6, col("ventricle.P")[f]/1e3; color=c, lw=2, label=lab)
    od = loadold(sch, o, dt)
    if od !== nothing
        on, OM = od; ocol(x) = OM[:, findfirst(==(x), on)]; ot = ocol("time"); of = (ot .>= 1.6) .& (mod.(1:length(ot), 20) .== 0)
        scatter!(pPo, ot[of], ocol("Po")[of]; marker=:circle, ms=4, markerstrokecolor=:black, markercolor=:white, label="")
        scatter!(pPi, ot[of], ocol("Pi")[of]; marker=:circle, ms=4, markerstrokecolor=:black, markercolor=:white, label="")
        scatter!(pPV, ocol("Vv")[of]*1e6, ocol("Pv")[of]/1e3; marker=:circle, ms=4, markerstrokecolor=:black, markercolor=:white, label="")
    end
end
plot!(pPV, mcol("Vv")[fin_m]*1e6, [(1) for _ in 1:count(fin_m)]*NaN; label="")
fig = plot(pPo, pPi, pQ, pPV; layout=(2, 2), size=(1300, 900), margin=5Plots.mm,
           plot_title="Closed-loop pipe, final cardiac cycle (hollow circles: original runs)")
savefig(fig, joinpath(out, "pipeClosedLoop_finalCycle.png"))
println("wrote ", out)
