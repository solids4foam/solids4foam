# Convergence and waveform plots for the coupled pipe-Windkessel ladder.
#   julia plotPipeWindkessel.jl [runs] [outdir]
# New results are lines with markers; the stored results of the original
# hand-written Windkessel couplings (verification/reference/pipeWindkessel) are hollow black circles.
ENV["GKSwstype"] = "100"
using Plots, Printf

runs = length(ARGS) >= 1 ? ARGS[1] : "runs"
out = length(ARGS) >= 2 ? ARGS[2] : "plots"
mkpath(out)

const R, RCH, C, A, W = 10.0, 1.0, 0.01, 1.0, 2pi
function exact(t)
    k = W*R*C
    pc = A*R/(1 + k^2)*(sin(W*t) - k*cos(W*t) + k*exp(-t/(R*C)))
    return pc, pc + RCH*A*sin(W*t)
end

function readcoupling(path)
    names = String[]; rows = Vector{Vector{Float64}}()
    for l in eachline(path)
        if startswith(l, "#")
            p = split(l[2:end]); if !isempty(p) && p[1] == "time"; names = String.(p); end
        elseif !isempty(strip(l)); push!(rows, parse.(Float64, split(l))); end
    end
    ix = findfirst(==("outlet.Ppatch"), names)
    return [r[1] for r in rows], [r[ix] for r in rows]
end

rmse(t, po) = sqrt(sum((po[i] - exact(t[i])[2])^2 for i in eachindex(t) if 0 < t[i] <= 1 + 1e-9)/count(x -> 0 < x <= 1 + 1e-9, t))

# --- collect
res = Dict{Tuple{String,Int},Vector{Tuple{Float64,Float64}}}()
wave = Dict{Tuple{String,Int,String},Tuple{Vector{Float64},Vector{Float64}}}()
for scheme in readdir(runs), ord in readdir(joinpath(runs, scheme)), dt in readdir(joinpath(runs, scheme, ord))
    f = joinpath(runs, scheme, ord, dt, "postProcessing", "circulationCoupling", "0", "coupling.dat")
    isfile(f) || continue
    t, po = readcoupling(f)
    push!(get!(res, (scheme, parse(Int, ord)), []), (parse(Float64, dt), rmse(t, po)))
    wave[(scheme, parse(Int, ord), dt)] = (t, po)
end

# --- stored results
old = Dict{Tuple{String,Int},Vector{Tuple{Float64,Float64}}}()
base = joinpath(@__DIR__, "..", "..", "verification", "reference", "pipeWindkessel")
ordname = Dict("First" => 1, "Second" => 2, "Third" => 3, "Fourth" => 4)
for f in ("orderAccuracy_first_second.csv", "orderAccuracy_third_fourth.csv")
    p = joinpath(base, f); isfile(p) || continue
    ls = readlines(p); h = split(ls[1], ",")
    ia, id, ir = findfirst(==("algorithm"), h), findfirst(==("deltaT"), h), findfirst(==("rmse"), h)
    for l in ls[2:end]
        c = split(l, ","); m = match(r"^(Explicit|Implicit)(\w+?)Order$", c[ia])
        m === nothing && continue
        push!(get!(old, (lowercase(m[1]), ordname[m[2]]), []), (parse(Float64, c[id]), parse(Float64, c[ir])))
    end
end

col = Dict(1 => :royalblue, 2 => :darkorange, 3 => :seagreen, 4 => :crimson)
conv = plot(; xscale=:log10, yscale=:log10, xlabel="time step dt (s)", ylabel="RMSE of outlet pressure Po (Pa)",
            title="Coupled pipe + Windkessel: temporal convergence", legend=:bottomright,
            framestyle=:box, size=(900, 650), minorgrid=true)
for ((scheme, ord), pts) in sort(collect(res); by=x -> (x[1][1], x[1][2]))
    sort!(pts)
    plot!(conv, first.(pts), last.(pts); color=col[ord], marker=(scheme == "explicit" ? :square : :diamond),
          ls=(scheme == "explicit" ? :dash : :solid), lw=2, label="$scheme, order $ord")
    if haskey(old, (scheme, ord))
        o = filter(p -> any(isapprox(p[1], q[1]; rtol=1e-9) for q in pts), old[(scheme, ord)])
        scatter!(conv, first.(o), last.(o); marker=:circle, ms=6, markerstrokecolor=:black,
                 markercolor=:white, markerstrokewidth=1.5, label="")
    end
end
scatter!(conv, [NaN], [NaN]; marker=:circle, ms=8, markerstrokecolor=:black, markercolor=:white, label="original hand-written couplings")
for (ord, y0) in ((1, 2.0), (2, 0.5), (3, 0.05), (4, 0.01))
    plot!(conv, [0.005, 0.1], [y0*(0.005/0.1)^ord, y0]; color=:gray, ls=:dot, lw=1, label="")
    annotate!(conv, 0.115, y0, text("dt^$ord", 8, :gray, :left))
end
savefig(conv, joinpath(out, "pipeWindkessel_convergence.png"))

# --- waveform at one dt
dt = "0.05"
ts = range(0, 1; length=400)
wp = plot(ts, [exact(t)[2] for t in ts]; color=:black, lw=3, label="exact", xlabel="t (s)", ylabel="Po (Pa)",
          title="Outlet pressure, dt = $dt s", framestyle=:box, size=(900, 600), legend=:topright)
for ((scheme, ord, d), (t, po)) in sort(collect(wave); by=x -> (x[1][1], x[1][2]))
    d == dt || continue
    plot!(wp, t, po; color=col[ord], marker=(scheme == "explicit" ? :square : :diamond), ms=4,
          ls=(scheme == "explicit" ? :dash : :solid), label="$scheme $ord")
end
savefig(wp, joinpath(out, "pipeWindkessel_outletPressure_dt0p05.png"))

# --- pointwise error at one dt
ep = plot(; xlabel="t (s)", ylabel="Po - exact (Pa)", title="Pointwise error, dt = $dt s", framestyle=:box,
          size=(900, 600), legend=:bottomleft)
for ((scheme, ord, d), (t, po)) in sort(collect(wave); by=x -> (x[1][1], x[1][2]))
    d == dt || continue
    plot!(ep, t, [po[i] - exact(t[i])[2] for i in eachindex(t)]; color=col[ord],
          ls=(scheme == "explicit" ? :dash : :solid), lw=2, label="$scheme $ord")
end
savefig(ep, joinpath(out, "pipeWindkessel_error_dt0p05.png"))
println("wrote ", out)
