# ==============================================================================
# File:        CascAnalytics.jl
# Project:     A Criterion for Safe Overshoot in Coupled Tipping Systems
#
# Description:
#   Module containing functions used in main.jl.
#
# Author(s):   Sacha Sinet (s.a.m.sinet@uu.nl)
# Affiliation: Institute for Marine and Atmospheric research Utrecht (IMAU)
#
# License:     MIT License (see LICENSE file)
#
# ==============================================================================

module CascAnalytics
using DifferentialEquations
using BifurcationKit
using UnPack

#region General Functions
function c(state, pvec, t) #coupling functions used for each models
    explabel = pvec[1]

    if explabel == 1
        x, y = state
        explabel, a, b, μ0, r, distance, j, ϵ, γ = pvec

        der = (3 * x - x^3 + μ0 + r * t) # careful if you use a rate, then 2 is not really the position
        return der

    elseif explabel == 2
        x, y = state
        explabel, a, b, μ0, r, distance, j, ϵ, γ = pvec

        coup = 1/5 * (cos(π/14 - 2π/7*x))^2

        return coup

    elseif explabel == 3
        x, y, P, T = state
        explabel, α, μcar, tdiff, F0, Q0, r, rP, Pd, b, K, hP, rm, mA, hA, mf, hf, β, P0, T0, distance, j, ϵ, γ = pvec

        ΔQ = 1 + μcar * (x - y)^2 - Q0
        return ΔQ



    elseif explabel == 4
        V, x, y = state
        explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j, ϵ, γ = pvec

        flux = -(alead * V^3 + blead * V^2 + clead * V + (Vp - Vm)^3 / (2 * (Tm - Tp)) * (Tf + r .* t) + (Tp * Vm^2 * (Vm - 3Vp) - Tm * Vp^2 * (Vp - 3Vm)) / (2 * (Tm - Tp))) # careful if you use a rate, then Tp is not really the position
        return flux
    end
end

function stopsimulation!(integrator) #callback for stopping simulation when tipping
    terminate!(integrator)
end

function Map(ϵrange, γrange, tspan, prob, prob_func, tippingcallback) #create Map in (ϵ,γ)
    param_list = [(p1, p2) for p1 in ϵrange, p2 in γrange]
    param_list = vec(param_list)

    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func=(sol, i) -> (sol.t[end], false))
    sol = solve(ensemble_prob, Rodas5(), EnsembleSerial(); trajectories=length(param_list), callback=tippingcallback)
    result = CascAnalytics.result_from_sol(sol, tspan, ϵrange, γrange)

    return result
end

function result_from_sol(sol, tspan, p1range, p2range) #Get result from the Map
    result = NaN .* zeros(Float64, length(p1range), length(p2range))
    for i ∈ 1:length(sol)
        sol[i] >= tspan[end] ? tip = 0 : tip = 1
        result[1+mod(i - 1, length(p1range)), 1+div(i - 1, length(p1range))] = tip
    end
    return result
end

function δ1f(statemax, pnupl, tmax) #Jacobian functions used for each models
    @unpack explabel, = pnupl
    if explabel == 1 || explabel == 2
        x, y = statemax
        δ1f = hcat((3 - 3 * x^2))
    elseif explabel == 3
        x, y = statemax
        @unpack α, μcar, tdiff, F0, r = pnupl
        δ1f = [[-1 - α - ((x - y)^2) * μcar - 2x * (x - y) * μcar, 2x * (x - y) * μcar] [-2(x - y) * y * μcar, -1 - ((x - y)^2) * μcar + 2(x - y) * y * μcar]]

    elseif explabel == 4
        V, x, y = statemax
        @unpack explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j = pnupl
        δ1f = hcat((3 * alead * V^2 + 2 * blead * V + clead))
    end

    return δ1f
end

function δ11c(statemax, pnupl, tmax) #Hessian functions used for each models
    @unpack explabel, = pnupl
    if explabel == 1
        x, y = statemax
        δ11c = hcat(-(-6 * x))
    elseif explabel == 2
        x, y = statemax
        δ11c = hcat(8/245*π^2*cos(2/7*π*(2*x+3)))
    elseif explabel == 3
        x, y = statemax
        @unpack α, μcar, tdiff, F0, r = pnupl
        δ11c = [[2μcar, -2μcar] [-2μcar, 2μcar]]

    elseif explabel == 4
        V, x, y = statemax
        @unpack explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j = pnupl
        δ11c = hcat(-(6 * alead * V + 2 * blead))
    end

    return δ11c
end

function δ1c(statemax, pnupl, tmax) #Gradient functions used for each models
    @unpack explabel, = pnupl
    if explabel == 1
        x, y = statemax
        δ1c = [-(3 - 3 * x^2),]
    elseif explabel == 2
        x, y = statemax
        δ1c = [2/35*π*sin(2/7*π*(2*x+3)),]
    elseif explabel == 3
        x, y = statemax
        @unpack α, μcar, tdiff, F0, r = pnupl
        δ1c = [2(x - y) * μcar, -2(x - y) * μcar]

    elseif explabel == 4
        V, x, y = statemax
        @unpack explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j = pnupl
        δ1c = [-(3 * alead * V^2 + 2 * blead * V + clead),]
    end

    return δ1c
end

function f(statemax, pnupl, tmax) #R.H.S for each models
    @unpack explabel, = pnupl
    if explabel == 1 || explabel == 2
        x, y = statemax
        @unpack explabel, a, b, μ0, r, distance, j = pnupl
        return [(3 * x - x^3 + μ0 + r * tmax),]

    elseif explabel == 3
        x, y = statemax
        @unpack α, μcar, tdiff, F0, r = pnupl

        f1 = -α * (x - 1) - x * (1 + μcar * (x - y)^2)
        f2 = F0 + r * tmax - y * (1 + μcar * (x - y)^2) #Note that here, we take the system at the criticality, but normally there should be the rate

        return [f1, f2]

    elseif explabel == 4
        V, x, y = statemax
        @unpack explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j = pnupl
        return [alead * V^3 + blead * V^2 + clead * V + (Vp - Vm)^3 / (2 * (Tm - Tp)) * (Tf + r .* tmax) + (Tp * Vm^2 * (Vm - 3Vp) - Tm * Vp^2 * (Vp - 3Vm)) / (2 * (Tm - Tp)),]
    end
end

function δ2f(statemax, pnupl, tmax) #Derivative of the R.H.S. w.r.t. the slow time ϵt for each models
    @unpack explabel, = pnupl
    if explabel == 1
        @unpack explabel, a, b, μ0, r, distance, j = pnupl
        δ2f = [r,]
    elseif explabel == 2
        @unpack explabel, a, b, μ0, r, distance, j = pnupl
        δ2f = [r,]
    elseif explabel == 3
        x, y = statemax
        @unpack α, μcar, tdiff, F0, r = pnupl
        δ2f = [0., r]

    elseif explabel == 4
        V, x, y = statemax
        @unpack explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j = pnupl
        δ2f = [(Vp - Vm)^3 / (2 * (Tm - Tp)) * r,]
    end

    return δ2f
end

function δ12c(statemax, pnupl, tmax) #This is zero in each cases considered here (due to coupling expressions)
    @unpack explabel, = pnupl
    if explabel == 1
        δ12c = [0.,]
    elseif explabel == 2
        δ12c = [0.,]
    elseif explabel == 3
        δ12c = [0., 0.]

    elseif explabel == 4
        δ12c = [0.,]
    end

    return δ12c
end

function δ22c(statemax, pnupl, tmax) #This is zero in each cases considered here (due to coupling expressions)
    @unpack explabel, = pnupl
    if explabel == 1
        δ22c = 0.
    elseif explabel == 2
        δ22c = 0.
    elseif explabel == 3
        δ22c = 0.

    elseif explabel == 4
        δ22c = 0.
    end

    return δ22c
end

function S(statemax, pnupl, tmax) #Leading system coefficient S
    return (f(statemax, pnupl, tmax))' * δ11c(statemax, pnupl, tmax) * f(statemax, pnupl, tmax) + (δ1c(statemax, pnupl, tmax))' * δ1f(statemax, pnupl, tmax) * f(statemax, pnupl, tmax) + (δ1c(statemax, pnupl, tmax))' * δ2f(statemax, pnupl, tmax) + 2 * (δ12c(statemax, pnupl, tmax))' * f(statemax, pnupl, tmax) + δ22c(statemax, pnupl, tmax)
end

function criterium(ϵrange, γrange, parameters, statemax, F; tmax=0.0) #Criterium for safe overshoot
    criterium = NaN .* zeros(Float64, length(ϵrange), length(γrange))

    @unpack j, distance = parameters.nupl

    cmax = c(statemax, parameters.vec, tmax)

    for i ∈ 1:length(ϵrange)
        for k ∈ 1:length(γrange)
            criterium[i, k] = F * sqrt.(1 / (γrange[k] * ϵrange[i]^(2 + j) * abs(S(statemax, parameters.nupl, tmax)))) * (γrange[k] * ϵrange[i]^(j) * cmax - distance)
        end
    end

    return criterium
end

function criteriumcorr(ϵtrange, γtrange, parameters, F) #Criterium corrollary for safe overshoot
    criteriumcorr = NaN .* zeros(Float64, length(ϵtrange), length(γtrange))

    @unpack j, distance = parameters.nupl

    for i ∈ 1:length(ϵtrange)
        for k ∈ 1:length(γtrange)
            criteriumcorr[i, k] = F / 4 * sqrt.(1 / (6 * γtrange[k] * ϵtrange[i]^3)) * (γtrange[k] * ϵtrange[i] * 4 - parameters.nupl.distance)
        end
    end

    return criteriumcorr
end

#endregion

#region Functions for Generic system
struct parametersGeneric
    vec::Vector
    nupl::NamedTuple

    function parametersGeneric(;
        explabel=NaN,
        a=2,
        b=0.2,
        μ0=9 / 8,
        r=0.1,
    )
        explabel == 1 ? j = 1 : j = 0

        # Parameters ----------
        nupl = (explabel=explabel,
            a=a,
            b=b,
            μ0=μ0,
            r=r,
            distance=b,
            j=j,
        )

        vec = [nupl.explabel,
            nupl.a, nupl.b,
            nupl.μ0,
            nupl.r,
            nupl.distance,
            nupl.j]
        new(vec, nupl)
    end
end

function FGeneric!(dstate, state, pvec, t)
    x, y = state
    explabel, a, b, μ0, r, distance, j, ϵ, γ = pvec

    dstate[1] = 3 * x - x^3 + μ0 + r * t
    dstate[2] = (a * y^2 - b + γ * ϵ^j * CascAnalytics.c(state, pvec, t)) * 1 / ϵ

    return nothing
end

function tipevent(u, t, integrator)
    return (u[2] > 0.5)
end

function affecttip!(integrator)
    terminate!(integrator)
end

tippingcallback = DiscreteCallback(tipevent, affecttip!; save_positions=(false, false))
dstate = SavedValues(Float64, Vector{Float64})
cb1 = SavingCallback((u, t, integrator) -> integrator(t, Val{1}), dstate)
cbs = CallbackSet(cb1, tippingcallback)
#endregion

#region Functions for CessiVeg
struct parametersCessiVeg
    vec::Vector
    nupl::NamedTuple

    function parametersCessiVeg(; #default value of all parameters
        explabel=3,
        α=3.6 * 10^3, #all that has to do with Cessi
        μcar=6.2,
        tdiff=180,
        F0=1.1,
        Q0=NaN,
        r=1e-4,
        rP=1.0, #all that has to do with Veg
        Pd=3.9,
        b=1.32,
        K=90.,
        hP=0.5,
        rm=0.3,
        mA=0.15,
        hA=10.,
        mf=0.11,
        hf=64.,
        β=7,
        P0=3.,
        T0=100.,
        distance=NaN,
        j=0)

        # Put in nupl ----------
        nupl = (explabel=explabel,
            α=α,
            μcar=μcar,
            tdiff=tdiff,
            F0=F0,
            Q0=Q0,
            r=r,
            rP=rP,
            Pd=Pd,
            b=b,
            K=K,
            hP=hP,
            rm=rm,
            mA=mA,
            hA=hA,
            mf=mf,
            hf=hf,
            β=β,
            P0=P0,
            T0=T0,
            distance=distance,
            j=j)

        # Put in vec ----------
        vec = [explabel, α, μcar, tdiff, F0, Q0, r, rP, Pd, b, K, hP, rm, mA, hA, mf, hf, β, P0, T0, distance, j]



        new(vec, nupl)
    end
end

function FCessiVeg!(dstate, state, pvec, t)
    x, y, P, T = state
    explabel, α, μcar, tdiff, F0, Q0, r, rP, Pd, b, K, hP, rm, mA, hA, mf, hf, β, P0, T0, distance, j, ϵ, γ = pvec

    dstate[1] = -α * (x - 1) - x * (1 + μcar * (x - y)^2)
    dstate[2] = F0 + r * t - y * (1 + μcar * (x - y)^2)

    ΔQ = c(state, pvec, t)

    dstate[3] = (rP * (Pd + b * T * T0 / K + γ * ϵ^j * ΔQ - P * P0)) / P0 * 1 / ϵ
    dstate[4] = (P * P0 / (hP + P * P0) * rm * T * T0 * (1 - T * T0 / K) - mA * T * T0 * hA / (T * T0 + hA) - mf * T * T0 * hf^β / (hf^β + (T * T0)^β)) / T0 * 1 / ϵ

    return nothing
end

function tipeventCessiVeg(u, t, integrator)
    return (u[4] < 0.6)
end
tippingcallbackCessiVeg = DiscreteCallback(tipeventCessiVeg, stopsimulation!; save_positions=(false, false))
#endregion

#region Functions for GISCessi
struct parametersGISCessi
    vec::Vector
    nupl::NamedTuple

    function parametersGISCessi(; #default value of all parameters
        explabel=4,
        α=3.6 * 10^3, #all that has to do with Cessi
        μcar=6.2,
        tdiff=180,
        F0=1.1,
        Q0=NaN,
        Tf=0.,
        Tp=1.52,
        Tm=0.3,
        Vp=0.77,
        Vm=0.3526554620064224,
        tmelt=470.0,
        alead=-1,
        blead=3 * (Vm + Vp) / 2,
        clead=-3 * Vm * Vp,
        dlead=(Vp - Vm)^3 / (2 * (Tm - Tp)) * Tf + (Tp * Vm^2 * (Vm - 3Vp) - Tm * Vp^2 * (Vp - 3Vm)) / (2 * (Tm - Tp)),
        r=1e-4,
        distance=NaN,
        j=1)

        # Put in nupl ----------
        nupl = (explabel=explabel,
            α=α,
            μcar=μcar,
            tdiff=tdiff,
            F0=F0,
            Q0=Q0,
            Tf=Tf,
            Tp=Tp,
            Tm=Tm,
            Vp=Vp,
            Vm=Vm,
            tmelt=tmelt,
            alead=-1,
            blead=3 * (Vm + Vp) / 2,
            clead=-3 * Vm * Vp,
            dlead=(Vp - Vm)^3 / (2 * (Tm - Tp)) * Tf + (Tp * Vm^2 * (Vm - 3Vp) - Tm * Vp^2 * (Vp - 3Vm)) / (2 * (Tm - Tp)),
            r=r,
            distance=distance,
            j=j)

        # Put in vec ----------
        vec = [explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead, r, distance, j]



        new(vec, nupl)
    end
end

function FGISCessi!(dstate, state, pvec, t)
    V, x, y = state
    explabel, α, μcar, tdiff, F0, Q0, Tf, Tp, Tm, Vp, Vm, tmelt, alead, blead, clead, dlead,  r, distance, j, ϵ, γ = pvec

    dstate[1] = (alead * V^3 + blead * V^2 + clead * V + (Vp - Vm)^3 / (2 * (Tm - Tp)) * (Tf + r * t) + (Tp * Vm^2 * (Vm - 3Vp) - Tm * Vp^2 * (Vp - 3Vm)) / (2 * (Tm - Tp)))

    flux = CascAnalytics.c(state, pvec, t)

    dstate[2] = (-α * (x - 1) - x * (1 + μcar * (x - y)^2)) * 1 / ϵ
    dstate[3] = (F0 + γ * ϵ^j * flux - y * (1 + μcar * (x - y)^2)) * 1 / ϵ

    return nothing
end

function tipeventGISCessi(u, t, integrator)
    return (u[3] > 1.0)
end

tippingcallbackGISCessi = DiscreteCallback(tipeventGISCessi, CascAnalytics.stopsimulation!; save_positions=(false, false))
#endregion

end