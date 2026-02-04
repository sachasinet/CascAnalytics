# ==============================================================================
# File:        main.jl
# Project:     A Criterion for Safe Overshoot in Coupled Tipping Systems
#
# Description:
#   Numerical methods used to produce the data of the article.
#   Note that, for numerical efficiency, all computations are performed 
#   in terms of the slow time ϵt.
#
# Author(s):   Sacha Sinet (s.a.m.sinet@uu.nl)
# Affiliation: Institution for Marine and Atmospheric research Utrecht (IMAU)
#
# Created:     YYYY-MM-DD
# Last updated: YYYY-MM-DD
#
# License:     MIT License (see LICENSE file)
#
# Reference:
#   Author et al. (YEAR), Title, Journal, DOI
#   OR arXiv:XXXX.XXXXX
#
# ==============================================================================

#region Starting
#Select environment
cd("./Documents/Projects/CascAnalytics")
using Pkg
Pkg.activate(".")
include("./CascAnalytics.jl")

#Add Packages and Module
using DifferentialEquations, DiffEqCallbacks
using CairoMakie
CairoMakie.activate!(inline=true)
using LaTeXStrings
using BifurcationKit, UnPack
using JLD2

#Load bifurcation diagrams
@load "./Data/brScurve.jld2"
@load "./Data/brsn.jld2"
@load "./Data/brCessi.jld2"
@load "./Data/brVeg.jld2"
@load "./Data/brGIS.jld2"

#define colors
darkgreen = colorant"#117733";
lightblue = colorant"#6699CC";
lightyellow = colorant"#EECC66";
darkred = colorant"#994455";
lightred = colorant"#EE99AA";
darkgreen = colorant"#117733";
darkblue = colorant"#004488";
darkyellow = colorant"#997700";
lightgreen = cgrad([darkgreen, :white], alpha=0.8)[0.4];
verylightblue = cgrad([lightblue, :white], alpha=0.8)[0.4];

palette = cgrad([lightred, lightyellow], 2, categorical=true);
#endregion

#Note: for numerical purposes, we express systems in terms of the slow time ϵt. All equations are modified accordingly.

#region Generic derivative
#We run one trajectory, to compute tmax and statemax (in terms of the slow time ϵt)
r = 0.1
parametersDer = CascAnalytics.parametersGeneric(explabel=1, r=r, μ0=9 / 8)

ϵ = 0.05
γ = 1.4
append!(parametersDer.vec, ϵ, γ)

state0 = [1.5, sqrt(parametersDer.nupl.b / parametersDer.nupl.a)]
tspan = (0, 4 / r)

prob = ODEProblem(CascAnalytics.FGeneric!, state0, tspan, parametersDer.vec)
@time trajectory = solve(prob, RK4(); callback=CascAnalytics.cbs);

t = trajectory.t
x = Matrix(trajectory)[1, :]
y = Matrix(trajectory)[2, :]
coupling = @. -(3 * x - x^3 - parametersDer.nupl.μ0 .- parametersDer.nupl.r .* t)

begin
    f = Figure()
    ax1 = Axis(f[1, 1], limits=(nothing, nothing, -3., 3.))
    lines!(ax1, brScurve.branch.param, brScurve.branch.x; color=:grey)
    lines!(ax1, parametersDer.nupl.μ0 .+ parametersDer.nupl.r .* t, x, color=:red)
    ax1.title = L"\text{Leading system trajectory}"
    ax1.xlabel = L"\mu"
    ax1.ylabel = L"x"

    ax2 = Axis(f[1, 2], limits=(0, nothing, -1.0, 1.0))
    lines!(ax2, t, ϵ .* coupling; color=:red)
    ax2.title = L"\text{Forcing}"
    ax2.xlabel = L"\text{time}"
    ax2.ylabel = L"\dot{x}"

    Label(f[2, 1:2], L"\text{Following system trajectory for different coupling}")

    ax3 = Axis(f[3, 1], limits=(0.0, 0.3, -0.5, 0.5))
    lines!(ax3, brsn.branch.param, brsn.branch.x; color=:grey)
    la = lines!(ax3, γ .* coupling .* ϵ, y)
    ax3.xlabel = L"\gamma \dot{x}"
    ax3.ylabel = L"y"

    ax4 = Axis(f[3, 2], limits=(0, nothing, -0.5, 0.5))
    lines!(ax4, t, y)
    ax4.xlabel = L"\text{time}"
    ax4.ylabel = L"y"
end
f

# @save "./Data/trajDer1p4.jld2" t x y coupling parametersDer ϵ γ

idmax = findmax(coupling)[2]
statemax = [x[idmax], y[idmax]] #system at tmax
tmax = t[idmax] #tmax

#Map
ϵrangeDer = range(0.001, 0.1, 400)
γrangeDer = range(0.001, 3.0, 400)

param_list = [(p1, p2) for p1 in ϵrangeDer, p2 in γrangeDer]
param_list = vec(param_list)

function prob_func(prob, i, repeat)
    par = CascAnalytics.parametersGeneric(explabel=1, r=r, μ0=9 / 8)
    ϵ = param_list[i][1]
    γ = param_list[i][2]
    append!(par.vec, ϵ, γ)

    remake(prob, p=par.vec)
end

resultDer = CascAnalytics.Map(ϵrangeDer, γrangeDer, tspan, prob, prob_func, CascAnalytics.cbs);

contourf(ϵrangeDer, γrangeDer, resultDer; levels=0.0:1.0:2.0)

#we compute the value of the coefficient F. Here, can be done analytically
a0 = 1.0
κ = 2.0
F = a0 * sqrt(2 * κ)

eig1 = zeros(Float64, length(brsn.branch.param))
d = zeros(Float64, length(brsn.branch.param))
for i ∈ 1:(brsn.specialpoint[2].idx-1)
    eig1[i] = real(brsn.eig[i].eigenvals[1])
    d[i] = (brsn.eig[i].eigenvals[1] .* brsn.eig[i].eigenvals[1]) / (brsn.specialpoint[2].param - brsn.branch.param[i])
end
coef = sqrt(d[brsn.specialpoint[1].idx-1] / 2) #Why not a nimus here?


#we compute the criterium
critDer = CascAnalytics.criterium(ϵrangeDer, γrangeDer, parametersDer, statemax, F; tmax=tmax)

#we compute the corrollary
ϵtrange = ϵrangeDer
γtrange = γrangeDer

critcorr = CascAnalytics.criteriumcorr(ϵtrange, γtrange, parametersDer, F)

#we plot the results figure
begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    contourf!(ax, ϵrangeDer, γrangeDer, resultDer; colormap=palette)
    contour!(ax, ϵrangeDer, γrangeDer, critDer; levels=[1,], color=:black, linewidth=3) #why the sqrt(2)?
    contour!(ax, ϵrangeDer, γrangeDer, critcorr; levels=[1,], color=:black, linewidth=3, linestyle=:dash) #why the sqrt(2)?
    # lines!(ax, ϵrange, 0.2 ./ (4 .* ϵrange), color=:green, linewidth=2)

    lines!(ax, ϵrangeDer, parametersDer.nupl.distance ./ (4 .* ϵrangeDer), linewidth=3, color=:black, linestyle=:dot)

    band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), parametersDer.nupl.distance ./ (CascAnalytics.c(statemax, parametersDer.vec, tmax) .* ϵrangeDer), color=:white)
    band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), parametersDer.nupl.distance ./ (CascAnalytics.c(statemax, parametersDer.vec, tmax) .* ϵrangeDer), color=lightblue)

    ylims!(ax, 0, 3)

    ax.xlabel = L"Timescale separation $\epsilon$"
    ax.ylabel = L"Coupling strength $\gamma$"
end
fig

# @save "./Data/trajDer1p4.jld2" t x y coupling parametersDer
# @save "./Data/ResultDerGen.jld2" resultDer ϵrangeDer γrangeDer critDer parametersDer r criteriumcorr xmax tmax
# @save "./Data/ResultDerGenr1p0.jld2" resultDer ϵrangeDer γrangeDer critDer parametersDer r criteriumcorr xmax tmax
#endregion

#region Generic nonlinear
#We run one trajectory, to compute tmax and statemax (in terms of the slow time ϵt)
r = 0.1
parametersNonl = CascAnalytics.parametersGeneric(explabel=2, r=r, μ0=9 / 8)

ϵ = 0.05
γ = 0.4
append!(parametersNonl.vec, ϵ, γ)

state0 = [1.5, sqrt(parametersNonl.nupl.b / parametersNonl.nupl.a)]
tspan = (0, 4 / r)

prob = ODEProblem(CascAnalytics.FGeneric!, state0, tspan, parametersNonl.vec)
@time trajectory = solve(prob, RK4(); callback=CascAnalytics.cbs);

t = trajectory.t
x = Matrix(trajectory)[1, :]
y = Matrix(trajectory)[2, :]
coupling = @. 0.2 * (cos((2 * π) / 7 * x + (15 * π) / 14))^2

begin
    f = Figure()
    ax1 = Axis(f[1, 1], limits=(nothing, nothing, -3., 3.))
    lines!(ax1, brScurve.branch.param, brScurve.branch.x; color=:grey)
    lines!(ax1, parametersNonl.nupl.μ0 .+ parametersNonl.nupl.r .* t, x, color=:red)
    ax1.title = L"\text{Leading system trajectory}"
    ax1.xlabel = L"\mu"
    ax1.ylabel = L"x"

    ax2 = Axis(f[1, 2], limits=(0, nothing, -1.0, 1.0))
    lines!(ax2, t, coupling; color=:red)
    ax2.title = L"\text{Forcing}"
    ax2.xlabel = L"\text{time}"
    ax2.ylabel = L"\dot{x}"

    Label(f[2, 1:2], L"\text{Following system trajectory for different coupling}")

    ax3 = Axis(f[3, 1], limits=(0.0, 0.3, -0.5, 0.5))
    lines!(ax3, brsn.branch.param, brsn.branch.x; color=:grey)
    la = lines!(ax3, γ .* coupling, y)
    ax3.xlabel = L"\gamma c_1(x)"
    ax3.ylabel = L"y"

    ax4 = Axis(f[3, 2], limits=(0, nothing, -0.5, 0.5))
    lines!(ax4, t, y)
    ax4.xlabel = L"\text{time}"
    ax4.ylabel = L"y"
end
f

idmax = findmax(coupling)[2]
statemax = [x[idmax], y[idmax]] #system at tmax
tmax = t[idmax] #tmax

#Map
ϵrangeNonl = range(0.001, 0.1, 400)
γrangeNonl = range(0.001, 3.0, 400)

param_list = [(p1, p2) for p1 in ϵrangeNonl, p2 in γrangeNonl]
param_list = vec(param_list)

function prob_func(prob, i, repeat)
    par = CascAnalytics.parametersGeneric(explabel=2, r=r, μ0=9 / 8)
    ϵ = param_list[i][1]
    γ = param_list[i][2]
    append!(par.vec, ϵ, γ)

    remake(prob, p=par.vec)
end

resultNonl = CascAnalytics.Map(ϵrangeNonl, γrangeNonl, tspan, prob, prob_func, CascAnalytics.cbs);

contourf(ϵrangeNonl, γrangeNonl, resultNonl; levels=0.0:1.0:2.0)

#we compute the value of the coefficient F. Here, can be done analytically
a0 = 1.0
κ = 2.0
F = a0 * sqrt(2 * κ)

#we compute the criterium
critNonl = CascAnalytics.criterium(ϵrangeNonl, γrangeNonl, parametersNonl, statemax, F; tmax=tmax)

#we plot the results figure
begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    contourf!(ax, ϵrangeNonl, γrangeNonl, resultNonl; colormap=palette)
    contour!(ax, ϵrangeNonl, γrangeNonl, critNonl; levels=[1,], color=:black, linewidth=3) #why the sqrt(2)?

    ylims!(ax, 0, 2)

    ax.xlabel = L"Timescale separation $\epsilon$"
    ax.ylabel = L"Coupling strength $\gamma$"
end
fig
#endregion

#region CessiVeg
#We run one trajectory, to compute tmax and statemax (in terms of the slow time ϵt)
r = 1e-4;
parametersCessiVeg = CascAnalytics.parametersCessiVeg(r=r, Q0=brCessi.branch.Q[1], distance=-(brVeg.branch.param[1] - brVeg.specialpoint[1].param));

ϵ = 1 / parametersCessiVeg.nupl.tdiff; #there should be the reference one
γ = 0.62;
append!(parametersCessiVeg.vec, ϵ, γ);

state0 = [0.9987362171488928, 0.24147297094357806, 5.0038820248103635 / parametersCessiVeg.nupl.P0, 75.5823037215545 / parametersCessiVeg.nupl.T0]; #equilibrium at reference
tspan = (0, (2.0 - 1.1) / r);

prob = ODEProblem(CascAnalytics.FCessiVeg!, state0, tspan, parametersCessiVeg.vec);
@time trajectory = solve(prob, RK4(), dt=0.004; abstol=1e-10, reltol=1e-8, maxiters=1e7);

t = trajectory.t;
x = Matrix(trajectory)[1, :];
y = Matrix(trajectory)[2, :];
P = Matrix(trajectory)[3, :];
T = Matrix(trajectory)[4, :];
Q = 1 .+ parametersCessiVeg.nupl.μcar .* (x .- y) .^ 2;

fig = Figure()
begin
    ax1 = Axis(fig[1, 1])
    lines!(ax1, brCessi.branch.param, brCessi.branch.Q)
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q)
    # scatter!(ax1, [parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t[overthresholdid]], [Q[overthresholdid]]; color=:red)
    # scatter!(ax1, [parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t[end - undertresholdid + 1]], [Q[end - undertresholdid + 1]]; color=:red)

    ax2 = Axis(fig[1, 2])
    lines!(ax2, brVeg.branch.param, brVeg.branch.T)
    lines!(ax2, parametersCessiVeg.nupl.Pd .+ γ .* (Q .- brCessi.branch.Q[1]), T)

    ax3 = Axis(fig[2, 1])
    lines!(ax3, t, Q)
    lines!(ax3, t, T)
end
fig

idmax = findmin(Q .- parametersCessiVeg.nupl.Q0)[2]
statemax = [x[idmax], y[idmax], NaN, NaN] #system at tmax
tmax = t[idmax] #tmax

#we find tiptimeCessi (for the reference value of ϵ and low r)
fluxx = -parametersCessiVeg.nupl.α .* (x .- 1) - x .* (1 .+ parametersCessiVeg.nupl.μcar .* (x .- y).^2);
fluxy = parametersCessiVeg.nupl.F0 .+ r .* t .- y .* (1 .+ parametersCessiVeg.nupl.μcar .* (x .- y).^2);
flux = sqrt.(fluxx .^ 2 .+ fluxy .^ 2);

threshold = 0.3;
overthresholdid = findfirst(x -> x > threshold, flux./maximum(flux));
undertresholdid = findfirst(x -> x > threshold, reverse(flux./maximum(flux)));
tiptimeCessi = (t[end - undertresholdid + 1] - t[overthresholdid])*180

#Map
ϵrangeCessiVeg = range(0.001, 0.045, 400) # normally upper bound at 1 ./50 * 1/180 * tiptimeAMOC = 0.044439423848619986
γrangeCessiVeg = range(0.0, 0.88, 400)

param_list = [(p1, p2) for p1 in ϵrangeCessiVeg, p2 in γrangeCessiVeg]
param_list = vec(param_list)

function prob_func(prob, i, repeat)
    par = CascAnalytics.parametersCessiVeg(r=r, Q0=brCessi.branch.Q[1], distance=-(brVeg.branch.param[1] - brVeg.specialpoint[1].param))
    ϵ = param_list[i][1]
    γ = param_list[i][2]
    append!(par.vec, ϵ, γ)

    remake(prob, p=par.vec)
end

resultCessiVeg = CascAnalytics.Map(ϵrangeCessiVeg, γrangeCessiVeg, tspan, prob, prob_func, CascAnalytics.tippingcallbackCessiVeg)

contourf(ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; levels=0.0:1.0:2.0)

#we compute the value of the coefficient F. Here, has to be computed numerically
eig1 = zeros(Float64, length(brVeg.branch.param))
eig2 = zeros(Float64, length(brVeg.branch.param))
d = zeros(Float64, length(brVeg.branch.param))
for i ∈ 1:(brVeg.specialpoint[1].idx-1)
    eig1[i] = real(brVeg.eig[i].eigenvals[1])
    eig2[i] = real(brVeg.eig[i].eigenvals[2])
    d[i] = (brVeg.eig[i].eigenvals[1] .* brVeg.eig[i].eigenvals[1]) / (brVeg.specialpoint[1].param - brVeg.branch.param[i])
end
F = sqrt(-d[brVeg.specialpoint[1].idx-1] / 2)

#we compute the criterium
critCessiVeg = CascAnalytics.criterium(ϵrangeCessiVeg, γrangeCessiVeg, parametersCessiVeg, statemax, F; tmax=tmax)

begin
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    lines!(ax1, brCessi.branch.param, brCessi.branch.Q)
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q)
    # scatter!(ax1, [parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t[overthresholdid]], [Q[overthresholdid]]; color=:red)
    # scatter!(ax1, [parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t[end - undertresholdid + 1]], [Q[end - undertresholdid + 1]]; color=:red)

    ax2 = Axis(fig[1, 2])
    lines!(ax2, brVeg.branch.param, brVeg.branch.T)
    lines!(ax2, parametersCessiVeg.nupl.Pd .+ γ .* (Q .- brCessi.branch.Q[1]), T)

    ax3 = Axis(fig[2, 1])
    lines!(ax3, t, Q)
    lines!(ax3, t, T)

    ax4 = Axis(fig[2:3, 2])
    heatmap!(ax4, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap=Reverse(:greys))
    contour!(ax4, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels=[-1,], color=:red, linewidth=1)
    hlines!([0.57])

    scatter!(ax4, ϵ, γ)

    ax4.xlabel = L"Timescale separation $\epsilon$"
    ax4.ylabel = L"Coupling strength $\gamma$"
end
fig
#endregion

#region GISCessi 
#We run one trajectory, to compute tmax and statemax (in terms of the slow time ϵt)
r = 0.5e-2;
parametersGISCessi = CascAnalytics.parametersGISCessi(r=r, Q0 = brCessi.branch.Q[1], Tf = 1.1, distance=brCessi.specialpoint[1].param - brCessi.param[1]);

ϵ = parametersGISCessi.nupl.tdiff/parametersGISCessi.nupl.τmelt;
γ = 3.8;
append!(parametersGISCessi.vec, ϵ, γ)

state0 = [0.8987549713444614, 0.9987362171488928, 0.24147297094357806];
tspan = (0, 2*0.9*1/r);

prob = ODEProblem(CascAnalytics.FGISCessi!, state0, tspan, parametersGISCessi.vec)
@time trajectory = solve(prob, RK4(), dt=0.004; abstol=1e-10, reltol=1e-8, maxiters=1e7);

t = trajectory.t;
V = Matrix(trajectory)[1, :];
x = Matrix(trajectory)[2, :];
y = Matrix(trajectory)[3, :];
Q = Q = 1 .+ parametersGISCessi.nupl.μcar .* (x .- y) .^ 2;
flux = @. -(parametersGISCessi.nupl.alead*V^3 + parametersGISCessi.nupl.blead*V^2 + parametersGISCessi.nupl.clead*V + (parametersGISCessi.nupl.Vp - parametersGISCessi.nupl.Vm)^3 / (2 * (parametersGISCessi.nupl.Tm - parametersGISCessi.nupl.Tp)) * (parametersGISCessi.nupl.Tf + parametersGISCessi.nupl.r * t) + (parametersGISCessi.nupl.Tp * parametersGISCessi.nupl.Vm^2 * (parametersGISCessi.nupl.Vm - 3parametersGISCessi.nupl.Vp) - parametersGISCessi.nupl.Tm * parametersGISCessi.nupl.Vp^2 * (parametersGISCessi.nupl.Vp - 3parametersGISCessi.nupl.Vm)) / (2 * (parametersGISCessi.nupl.Tm - parametersGISCessi.nupl.Tp)));

begin
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    lines!(ax1, brGIS.branch.param, brGIS.branch.V; color = :grey)
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :blue)

    ax2 = Axis(fig[1, 2])
    lines!(ax2, brCessi.branch.param, brCessi.branch.Q; color = :grey)
    lines!(ax2, parametersGISCessi.nupl.F0 .+ γ .* ϵ .* flux, Q; color = :green)

    ax3 = Axis(fig[2, 1])
    lines!(ax3, parametersGISCessi.nupl.τmelt .* t, V; color = :blue)
    lines!(ax3, parametersGISCessi.nupl.τmelt .* t, y; color = :green)

    xlims!(ax3,2000000, 2100000)
end
fig

# @save "./Data/trajDer1p4.jld2" t x y coupling parametersDer ϵ γ

idmax = findmax(flux)[2]
statemax = [V[idmax], NaN, NaN] #system at tmax
tmax = t[idmax] #tmax

#we find tiptimeGIS (for the reference value of ϵ and low r)
threshold = 0.3;
overthresholdid = findfirst(x -> x > threshold, flux./maximum(flux));
undertresholdid = findfirst(x -> x > threshold, reverse(flux./maximum(flux)));
tiptimeGIS = (t[end - undertresholdid + 1] - t[overthresholdid])*470

#Map
ϵrangeGISCessi = range(0.001, 3.1, 400) #normally 1 ./1000 * 180/470 * tiptimeGIS = 3.052716098349083
γrangeGISCessi = range(0.01, 15.0, 400)

param_list = [(p1, p2) for p1 in ϵrangeGISCessi, p2 in γrangeGISCessi]
param_list = vec(param_list)

function prob_func(prob, i, repeat)
    par = CascAnalytics.parametersGISCessi(r=r, Tf = 1.1)
    ϵ = param_list[i][1]
    γ = param_list[i][2]
    append!(par.vec, ϵ, γ)

    remake(prob, p=par.vec)
end

resultGISCessi = CascAnalytics.Map(ϵrangeGISCessi, γrangeGISCessi, tspan, prob, prob_func, CascAnalytics.tippingcallbackGISCessi);

contourf(ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; levels=0.0:1.0:2.0)

#we compute the value of the coefficient F. Here, computed numerically
eig1 = zeros(Float64, length(brCessi.branch.param));
eig2 = zeros(Float64, length(brCessi.branch.param));
d = zeros(Float64, length(brCessi.branch.param));
for i ∈ 1:(brCessi.specialpoint[1].idx-1)
    eig1[i] = real(brCessi.eig[i].eigenvals[1])
    eig2[i] = real(brCessi.eig[i].eigenvals[2])
    d[i] = (brCessi.eig[i].eigenvals[1] .* brCessi.eig[i].eigenvals[1]) / (brCessi.specialpoint[1].param - brCessi.branch.param[i])
end
F = sqrt(d[brCessi.specialpoint[1].idx-1] / 2)


#we compute the criterium
critGISCessi = CascAnalytics.criterium(ϵrangeGISCessi, γrangeGISCessi, parametersGISCessi, statemax, F; tmax=tmax)

#we compute the corrollary
a3 = parametersGISCessi.nupl.alead
a2 = parametersGISCessi.nupl.blead
a1 = parametersGISCessi.nupl.clead
a0 = parametersGISCessi.nupl.dlead #this in somputed at the critical, as has to be

ϵtrange = ϵrangeGISCessi .* (3*a3*a1-a2^2)/(9*a3)
γtrange = γrangeGISCessi .* -1/(3*a3)*(a2^2 - 3*a3*a1)^(1/2)

critcorr = CascAnalytics.criteriumcorr(ϵtrange, γtrange, parametersGISCessi, F)

#we plot the results figure
begin
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    lines!(ax1, brGIS.branch.param, brGIS.branch.V; color = :grey)
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :blue)

    ax2 = Axis(fig[1, 2])
    lines!(ax2, brCessi.branch.param, brCessi.branch.Q; color = :grey)
    lines!(ax2, parametersGISCessi.nupl.F0 .+ γ .* ϵ .* flux, Q; color = :green)

    ax3 = Axis(fig[2, 1])
    lines!(ax3, parametersGISCessi.nupl.τmelt .* t, V; color = :blue)
    lines!(ax3, parametersGISCessi.nupl.τmelt .* t, y; color = :green)

    xlims!(ax3,2000000, 2100000)

    ax4 = Axis(fig[2, 2])

    contourf!(ax4, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi)
    contour!(ax4, ϵrangeGISCessi, γrangeGISCessi, critGISCessi; levels=[1,], color=:black, linewidth=3, linestyle = :dash) #why the sqrt(2)?
    contour!(ax4, ϵrangeGISCessi, γrangeGISCessi, critcorr; levels=[1,], color=:black, linewidth=3, linestyle = :dash) #why the sqrt(2)?
    # lines!(ax, ϵrange, 0.2 ./ (4 .* ϵrange), color=:green, linewidth=2)

    band!(ax4, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi))
    band!(ax4, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi))

    scatter!(ax4, ϵ, γ)

    ylims!(ax4, 0 ,15)

    ax4.xlabel = L"Timescale separation $\epsilon$"
    ax4.ylabel = L"Coupling strength $\gamma$"
end
fig

#we compute the range of safe overshoot for the reference γ
γrangeGISCessi[101] #where the physical coupling is more or Less
couplingcomp = CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi .* γrangeGISCessi[101] #the coupling component for this γ
noid = findlast(x -> x < parametersGISCessi.nupl.distance, couplingcomp) #where no overshoot stops
safeid = findfirst(x -> x >= 1.0, resultGISCessi[:,101]) #where safe overshoot stops

ϵno =  ϵrangeGISCessi[noid]
ϵsafe =  ϵrangeGISCessi[safeid]

tiptimeGISno = 1 ./ ϵno * 180/470 * tiptimeGIS
tiptimeGISsafe = 1 ./ ϵsafe * 180/470 * tiptimeGIS

saferange = tiptimeGISno - tiptimeGISsafe
#endregion