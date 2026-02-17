# ==============================================================================
# File:        Plotting.jl
# Project:     A Criterion for Safe Overshoot in Coupled Tipping Systems
#
# Description:
#   Script producing all figures of the article.
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

#region Select environment
cd("./Documents/Projects/CascAnalytics")
using Pkg
Pkg.activate(".")
include("./CascAnalytics.jl")
#endregion

#region Add Packages and Module
using DifferentialEquations, DiffEqCallbacks
using CairoMakie
CairoMakie.activate!(inline=true)
using LaTeXStrings
using BifurcationKit, UnPack
using ColorSchemes
using JLD2
#endregion

#region constants
darkgreen = colorant"#117733";
lightblue = colorant"#0077BB";
lightyellow = colorant"#EECC66";
darkred = colorant"#994455";
lightred = colorant"#EE99AA";
darkgreen = colorant"#117733";
darkblue = colorant"#004488";
darkyellow = colorant"#997700";
lightgreen = cgrad([darkgreen, :white], alpha=0.8)[0.4];
verylightblue = cgrad([lightblue, :white], alpha=0.8)[0.4];

bluepalette = cgrad([:blue, :white], 100, categorical=true);

nocolor = ColorSchemes.seaborn_colorblind[10];
safecolor = ColorSchemes.seaborn_colorblind[3];
unsafecolor = ColorSchemes.seaborn_colorblind[2];
CessiVegcolor = ColorSchemes.seaborn_colorblind[4];
nonlcolor = ColorSchemes.seaborn_colorblind[1];
leadcolor = ColorSchemes.seaborn_colorblind[8];
#endregion

#region Fig1
# @load "./Data/cont.jld2"
@load "./Data/brScurve.jld2"
@load "./Data/brsn.jld2"

function addaxLabel(ax, label, offset; color=:black)
    xlow = ax.limits[][1]
    xup = ax.limits[][2]
    ylow = ax.limits[][3]
    yup = ax.limits[][4]

    posx = xlow + offset[1] * (xup - xlow)
    posy = yup - offset[2] * (yup - ylow)
    text!(ax, posx, posy; text=label, align=[:left, :center], fontsize=7, color=color, font="bold")
end

begin
    size_inches = (5.0, 5.0)
    size_pt = 72 .* size_inches

    f  = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Left column
    # -------------------------------------------------------------------------
    ga = f[1:2, 1] = GridLayout()
    Label(ga[1, 1, TopLeft()], L"(\text{a})";
        fontsize = 9,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    # --- Axes: trajectories ---------------------------------------------------
    @load "./Data/trajNonl0p8.jld2" t x y coupling parametersNonl ϵ

    ax1 = Axis(ga[1, 1], limits = (0, 20/ϵ, -2.3, 2.5))
    lines!(ax1, t/ϵ, x; color = leadcolor)
    ax1.title  = "Trajectories"
    ax1.ylabel = L"x"
    hidexdecorations!(ax1)
    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    addaxLabel(ax1, L"\text{Slow subsystem}", [0.02, 0.07])

    # --- Axes: coupling functions --------------------------------------------
    ax2a = Axis(ga[2, 1], limits = (0, 20/ϵ, -0.02, 0.3))

    # c1(x)
    @load "./Data/trajNonl0p8.jld2" t x y coupling parametersNonl ϵ
    lines!(ax2a, t/ϵ, coupling; color = nonlcolor)

    # c2(x)
    @load "./Data/trajDer0p8.jld2" t x y coupling parametersDer ϵ γ
    lines!(ax2a, t/ϵ, ϵ .* coupling; color = CessiVegcolor)

    hidexdecorations!(ax2a)
    hideydecorations!(ax2a, ticks = false, ticklabels = false, label = false)
    addaxLabel(ax2a, L"\text{Coupling functions}", [0.02, 0.07])
    addaxLabel(ax2a, L"c_2(x)", [0.67, 0.5]; color = CessiVegcolor)
    addaxLabel(ax2a, L"c_1(x)", [0.41, 0.5]; color = nonlcolor)

    # --- Axes: fast subsystem, y2 (via c2) -----------------------------------
    ax3 = Axis(ga[4, 1], limits = (0, 20/ϵ, -0.3, 0.6))

    @load "./Data/trajDer0p8.jld2" t x y coupling parametersDer ϵ γ
    lines!(ax3, t/ϵ, y; color = nocolor)

    @load "./Data/trajDer1p1.jld2" t x y coupling parametersDer ϵ γ
    lines!(ax3, t/ϵ, y; color = safecolor)

    @load "./Data/trajDer1p4.jld2" t x y coupling parametersDer γ
    lines!(ax3, t/ϵ, y; color = unsafecolor)

    ax3.ylabel = L"y_2"
    hidexdecorations!(ax3, ticks = false, ticklabels = false, label = false)
    hideydecorations!(ax3, ticks = false, ticklabels = false, label = false)
    addaxLabel(ax3, L"Coupling via $c_2(x)$", [0.02, 0.07])

    println(t[end])

    # --- Axes: fast subsystem, y1 (via c1) -----------------------------------
    ax4 = Axis(ga[3, 1], limits = (0, 20/ϵ, -0.3, 0.6))

    @load "./Data/trajNonl0p8.jld2" t x y coupling parametersNonl
    lines!(ax4, t/ϵ, y; color = nocolor)

    @load "./Data/trajNonl1p1.jld2" t x y coupling parametersNonl
    lines!(ax4, t/ϵ, y; color = safecolor)

    @load "./Data/trajNonl1p4.jld2" t x y coupling parametersNonl
    lines!(ax4, t/ϵ, y; color = unsafecolor)

    ax4.ylabel = L"y_1"
    ax4.xlabel = L"t"
    hideydecorations!(ax4, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax4)
    addaxLabel(ax4, L"Coupling via $c_1(x)$", [0.02, 0.07])

    println(t[end])

    rowgap!(ga, 0)

    # -------------------------------------------------------------------------
    # (b) Top-right: slow subsystem bifurcation + trajectory
    # -------------------------------------------------------------------------
    gb = f[1, 2] = GridLayout()
    Label(gb[1, 1, TopLeft()], L"(\text{b})";
        fontsize = 9,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    @load "./Data/trajNonl0p8.jld2" t x y coupling parametersNonl

    ax1 = Axis(
        gb[1, 1:2],
        limits = (
            parametersDer.nupl.μ0,
            3,
            nothing, nothing
        ),
        xticks = [0,1.5,2,2.5]
    )

    ibegin = brScurve.specialpoint[1].idx
    ibif1  = brScurve.specialpoint[2].idx
    ibif2  = brScurve.specialpoint[3].idx
    iend   = brScurve.specialpoint[4].idx

    lines!(ax1, brScurve.branch.param[ibegin:ibif1], brScurve.branch.x[ibegin:ibif1]; color = :black)
    lines!(ax1, brScurve.branch.param[ibif1:ibif2],  brScurve.branch.x[ibif1:ibif2];  color = :black, linestyle = :dash)
    lines!(ax1, brScurve.branch.param[ibif2:iend],   brScurve.branch.x[ibif2:iend];   color = :black)

    @load "./Data/trajDer0p8.jld2" t x y coupling parametersDer ϵ γ
    lines!(ax1, parametersDer.nupl.μ .+ parametersDer.nupl.r .* t, x; color = leadcolor)

    ax1.title  = "Slow subsystem"
    ax1.xlabel = L"\mu"
    ax1.ylabel = L"x"
    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    # -------------------------------------------------------------------------
    # (c) Bottom-right: fast subsystem phase portraits + legend
    # -------------------------------------------------------------------------
    gc = f[2, 2] = GridLayout()
    Label(gc[1, 1, TopLeft()], L"(\text{c})";
        fontsize = 9,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ibegin = brsn.specialpoint[1].idx
    ibif   = brsn.specialpoint[2].idx
    iend   = brsn.specialpoint[3].idx

    ax2 = Axis(gc[1, 1], limits = (0.0, 0.4, -0.4, 0.75), yticks = [0.0, 0.5], xticks = [0.1, 0.3])
    lines!(ax2, brsn.branch.param[ibegin:ibif], brsn.branch.x[ibegin:ibif]; color = :black)
    lines!(ax2, brsn.branch.param[ibif:iend],   brsn.branch.x[ibif:iend];   color = :black, linestyle = :dash)
    text!(ax2, -0.68, 1.85;
        text    = L"Coupling via $c_1(x)$",
        align   = [:left, :center],
        fontsize= 6,
        color   = :black,
        font    = "bold"
    )

    @load "./Data/trajNonl0p8.jld2" t x y coupling parametersNonl
    lines!(ax2, 0.8 .* coupling, y; color = nocolor)
    @load "./Data/trajNonl1p1.jld2" t x y coupling parametersNonl
    lines!(ax2, 1.2 .* coupling, y; color = safecolor)
    @load "./Data/trajNonl1p4.jld2" t x y coupling parametersNonl
    lines!(ax2, 1.4 .* coupling, y; color = unsafecolor)

    addaxLabel(ax2, L"Coupling via $c_1(x)$", [0.02, 0.05])
    ax2.xlabel = L"\gamma c_1(x)"
    ax2.ylabel = L"y_1"
    hideydecorations!(ax2, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax2, ticks = false, ticklabels = false, label = false)

    ax3 = Axis(gc[1, 2], limits = (0.0, 0.4, -0.4, 0.75), yticks = [0.0, 0.5], xticks = [0.1, 0.3])
    lines!(ax3, brsn.branch.param[ibegin:ibif], brsn.branch.x[ibegin:ibif]; color = :black)
    lines!(ax3, brsn.branch.param[ibif:iend],   brsn.branch.x[ibif:iend];   color = :black, linestyle = :dash)

    @load "./Data/trajDer0p8.jld2" t x y coupling parametersDer ϵ γ
    la = lines!(ax3, 0.8 * ϵ .* coupling, y; color = nocolor)
    @load "./Data/trajDer1p1.jld2" t x y coupling parametersDer ϵ γ
    lb = lines!(ax3, 1.1 * ϵ .* coupling, y; color = safecolor)
    @load "./Data/trajDer1p4.jld2" t x y coupling parametersDer ϵ γ
    lc = lines!(ax3, 1.4 * ϵ .* coupling, y; color = unsafecolor)

    addaxLabel(ax3, L"\text{Coupling via $c_2(x)$}", [0.02, 0.05])

    ax3.yaxisposition = :right
    ax3.xlabel = L"\gamma c_2(x)"
    ax3.ylabel = L"y_2"
    hideydecorations!(ax3, ticks = false, ticklabels = true, label = false)
    hidexdecorations!(ax3, ticks = false, ticklabels = false, label = false)

    ax4 = Axis(gc[1, 1:2], limits = (nothing, nothing, nothing, nothing))
    hidedecorations!(ax4)
    hidespines!(ax4)
    ax4.title = "Fast subsystem"

    Legend(
        f[3, 1:2],
        [la, lb, lc],
        [L"\gamma=0.8", L"\gamma=1.1", L"\gamma=1.4"];
        framevisible = false,
        orientation  = :horizontal,
        padding      = (0, 0, 0, 0)
    )

    colgap!(gc, 0)

    colgap!(f.layout, 1, Relative(0.04))
    rowgap!(f.layout, 1, Relative(0.03))
    rowgap!(f.layout, 2, Relative(0.03))

    f
end

save("./Manuscript/Figures/Fig1.pdf", f, pt_per_unit=1)
#endregion

#region Fig2

palette = cgrad([safecolor, unsafecolor], 2, categorical=true)

begin
    size_inches = (5.0, 3.0)
    size_pt     = 72 .* size_inches

    f  = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Nonlinear coupling (left panel)
    # -------------------------------------------------------------------------
    gb = f[1, 1] = GridLayout()
    Label(gb[1, 1, TopLeft()], L"(\text{a})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax = Axis(gb[1, 1])

    @load "./Data/ResultNonlGen.jld2" resultNonl ϵrangeNonl γrangeNonl critNonl parametersNonl r statemax tmax

    contourf!(ax, ϵrangeNonl, γrangeNonl, resultNonl; colormap = palette, levels = 2)
    contour!(ax, ϵrangeNonl, γrangeNonl, critNonl; levels = [1], color = :black, linewidth = 1)

    band!(ax, ϵrangeNonl, zeros(length(ϵrangeNonl)), 1.0; color = :white)
    band!(ax, ϵrangeNonl, zeros(length(ϵrangeNonl)), 1.0; color = nocolor)

    scatter!(ax, [0.05], [0.8]; marker = :star5, strokewidth = 1, markersize = 7, color = nocolor)
    scatter!(ax, [0.05], [1.1]; marker = :star5, strokewidth = 1, markersize = 7, color = safecolor)
    scatter!(ax, [0.05], [1.4]; marker = :star5, strokewidth = 1, markersize = 7, color = unsafecolor)

    ax.xlabel = L"Timescale separation $\epsilon$"
    ax.ylabel = L"Coupling strength $\gamma$"
    ax.title  = "Nonlinear coupling"

    # -------------------------------------------------------------------------
    # (b) Derivative coupling (right panel)
    # -------------------------------------------------------------------------
    ga = f[1, 2] = GridLayout()
    Label(ga[1, 1, TopLeft()], L"(\text{b})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax = Axis(ga[1, 1])

    @load "./Data/ResultDerGen.jld2" resultDer ϵrangeDer γrangeDer critDer parametersDer r critcorr statemax tmax

    contourf!(ax, ϵrangeDer, γrangeDer, resultDer; colormap = palette, levels = 2)
    contour!(ax, ϵrangeDer, γrangeDer, critDer;        levels = [1], color = :black, linewidth = 1)
    contour!(ax, ϵrangeDer, γrangeDer, critcorr;  levels = [1], color = :black, linewidth = 1, linestyle = :dash)

    band_y = parametersDer.nupl.distance ./ (CascAnalytics.c(statemax, parametersDer.vec, tmax) .* ϵrangeDer)
    band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), band_y; color = :white)
    band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), band_y; color = nocolor)

    scatter!(ax, [0.05], [0.8]; marker = :star5, strokewidth = 1, markersize = 7, color = nocolor)
    scatter!(ax, [0.05], [1.1]; marker = :star5, strokewidth = 1, markersize = 7, color = safecolor)
    scatter!(ax, [0.05], [1.4]; marker = :star5, strokewidth = 1, markersize = 7, color = unsafecolor)

    ylims!(ax, 0, 3)

    ax.xlabel = L"Timescale separation $\epsilon$"
    ax.ylabel = L"Coupling strength $\gamma$"
    ax.title  = "Derivative coupling"

    f
end

save("./Manuscript/Figures/Fig2.pdf", f, pt_per_unit=1)
#endregion

#region Fig3
palette = cgrad([safecolor, unsafecolor], 2, categorical=true)

begin
    size_inches = (5.0, 2.5)
    size_pt     = 72 .* size_inches

    f = Figure(size = size_pt, fontsize = 7, padding = 0)

    # Small helper to avoid repeating the same block 4×
    function add_panel!(ax, file; title, hidey = false)
        @load file resultDer ϵrangeDer γrangeDer critDer parametersDer r critcorr statemax tmax

        contourf!(ax, ϵrangeDer, γrangeDer, resultDer; colormap = palette, levels = 2)
        contour!(ax, ϵrangeDer, γrangeDer, critDer;       levels = [1], color = :black, linewidth = 1)
        contour!(ax, ϵrangeDer, γrangeDer, critcorr; levels = [1], color = :black, linewidth = 1, linestyle = :dash)

        band_y = parametersDer.nupl.distance ./ (CascAnalytics.c(statemax, parametersDer.vec, tmax) .* ϵrangeDer)
        band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), band_y; color = :white)
        band!(ax, ϵrangeDer, zeros(length(ϵrangeDer)), band_y; color = nocolor)

        ylims!(ax, 0, 3)

        ax.xlabel = L"Timescale separation $\epsilon$"
        ax.ylabel = L"Coupling strength $\gamma$"
        ax.title  = title

        if hidey
            hideydecorations!(ax)
        end

        return nothing
    end

    ax = Axis(f[1, 4], xticks = 0.02:0.04:0.08)
    add_panel!(ax, "./Data/ResultDerGenr1p0.jld2"; title = L"r=1", hidey = true)
    ax.ylabel = "" 

    ax = Axis(f[1, 3], xticks = 0.02:0.04:0.08)
    add_panel!(ax, "./Data/ResultDerGen.jld2"; title = L"r=10^{-1}", hidey = true)
    ax.ylabel = "" 

    ax = Axis(f[1, 2], xticks = 0.02:0.04:0.08)
    add_panel!(ax, "./Data/ResultDerGenr0p01.jld2"; title = L"r=10^{-2}", hidey = true)
    ax.ylabel = "" 

    ax = Axis(f[1, 1], xticks = 0.02:0.04:0.08)
    add_panel!(ax, "./Data/ResultDerGenr0p001.jld2"; title = L"r=10^{-3}", hidey = false)

    colgap!(f.layout, 0.0)

    f
end

save("./Manuscript/Figures/Fig3.pdf", f, pt_per_unit=1)
#endregion

#region Fig4
@load "./Data/brCessi.jld2"
@load "./Data/brVeg.jld2"

palette = cgrad([safecolor, unsafecolor], 2, categorical=true)

begin
    size_inches = (5.3, 3.7)
    size_pt     = 72 .* size_inches

    f  = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Time series (top-left)
    # -------------------------------------------------------------------------
    ga = f[1, 1] = GridLayout()
    Label(ga[1, 1, TopLeft()], L"(\text{a})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    axm1 = Axis(ga[1, 1:2], yticks = [1, 3, 5])

    @load "./Data/TrajCessiVegg0.62.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(axm1, t .* 180, Q; color = leadcolor)

    hidexdecorations!(axm1)
    hideydecorations!(axm1, ticks = false, ticklabels = false, label = false)
    ylims!(axm1, 0.7, 5.0)
    xlims!(axm1, 0, 600000)
    axm1.ylabel = L"Q"

    ax0 = Axis(ga[2, 1:2], yticks = [0.5, 0.7], xticks = [0, 200000, 400000, 600000])

    @load "./Data/TrajCessiVegg0.62.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax0, t .* 180, T; color = unsafecolor)

    @load "./Data/TrajCessiVegg0.5.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax0, t .* 180, T; color = nocolor)

    @load "./Data/TrajCessiVegg0.57.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax0, t .* 180, T; color = safecolor)

    hideydecorations!(ax0, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax0, ticks = false, ticklabels = false, label = false)
    ylims!(ax0, 0.4, 0.82)
    xlims!(ax0, 0, 600000)

    ax0.ylabel = L"T"
    ax0.xlabel = L"\text{Model years}"

    rowgap!(ga, 0)
    colgap!(ga, 0)

    # -------------------------------------------------------------------------
    # (b) Bifurcation-style panels (bottom-left)
    # -------------------------------------------------------------------------
    gb = f[2, 1] = GridLayout()
    Label(gb[1, 1, TopLeft()], L"(\text{b})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax1 = Axis(gb[1, 1]; xticks = [1.2, 1.4])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax1, brCessi.branch.param[1:sn1],   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax1, brCessi.branch.param[sn1:sn2], brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brCessi.branch.param[sn2:end], brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/TrajCessiVegg0.62.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q; color = leadcolor)

    xlims!(1.1, 1.4)

    ax1.ylabel = L"Q"
    ax1.xlabel = L"F"

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax2 = Axis(gb[1, 2], yaxisposition = :right, xticks = [-1.6, -0.8, 0.0])

    iend = brVeg.specialpoint[3].idx
    sn1  = brVeg.specialpoint[1].idx
    sn2  = brVeg.specialpoint[2].idx
    sn3  = brVeg.specialpoint[3].idx
    sn4  = brVeg.specialpoint[4].idx

    lines!(ax2, brVeg.branch.param[1:sn1]   .- parametersCessiVeg.nupl.Pd, brVeg.branch.T[1:sn1];   color = :black)
    lines!(ax2, brVeg.branch.param[sn1:sn2] .- parametersCessiVeg.nupl.Pd, brVeg.branch.T[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax2, brVeg.branch.param[sn2:sn3] .- parametersCessiVeg.nupl.Pd, brVeg.branch.T[sn2:sn3]; color = :black)

    @load "./Data/TrajCessiVegg0.62.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax2, γ .* (Q .- brCessi.branch.Q[1]), T; color = unsafecolor)

    @load "./Data/TrajCessiVegg0.57.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax2, γ .* (Q .- brCessi.branch.Q[1]), T; color = safecolor)

    @load "./Data/TrajCessiVegg0.5.jld2" r parametersCessiVeg ϵ γ state0 tspan prob t x y P T Q
    lines!(ax2, γ .* (Q .- brCessi.branch.Q[1]), T; color = nocolor)

    ax2.xlabel = L"\gamma\Delta Q"
    ax2.ylabel = L"T"

    xlims!(ax2, -2.3, -0.0)
    ylims!(ax2, 0.4, 0.8)

    hideydecorations!(ax2, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax2, ticks = false, ticklabels = false, label = false)

    colgap!(gb, 0)

    # -------------------------------------------------------------------------
    # (c) Parameter plane (right column)
    # -------------------------------------------------------------------------
    gc = f[1:2, 2] = GridLayout()
    Label(gc[1, 1, TopLeft()], L"(\text{c})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax3 = Axis(gc[1, 1])

    @load "./Data/ResultCessiVeg.jld2" resultCessiVeg ϵrangeCessiVeg γrangeCessiVeg critCessiVeg parametersCessiVeg tiptimeCessi statemax

    contourf!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap = palette, levels = 2)
    contour!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels = [-1], color = :black, linewidth = 1)

    γcrit = findfirst(x -> x > -1, critCessiVeg)

    band_y = parametersCessiVeg.nupl.distance / CascAnalytics.c(statemax, parametersCessiVeg.vec, 0)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = :white)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = nocolor)

    ylims!(ax3, 0.0, 0.88)

    vlines!(ax3, 0.005555; color = :grey, linestyle = :dot)

    scatter!(ax3, 0.005555, 0.62; marker = :star5, strokewidth = 1, markersize = 7, color = unsafecolor)
    scatter!(ax3, 0.005555, 0.57; marker = :star5, strokewidth = 1, markersize = 7, color = safecolor)
    scatter!(ax3, 0.005555, 0.5;  marker = :star5, strokewidth = 1, markersize = 7, color = nocolor)

    ax3.xlabel = L"Timescale separation $\epsilon$"
    ax3.ylabel = L"Coupling strength $\gamma$"

    # --- Top axis with AMOC tipping duration labels --------------------------
    ax3_top = Axis(gc[1, 1];
        xaxisposition = :top,
        yaxisposition = :right,
        xlabel        = L"\text{AMOC tipping duration (years)}",
    )

    # Make the top axis decorative only (keep your original behavior)
    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [400, 200, 100, 50]
    eps_ticks    = 1 ./ tiptimeticks * 1 / 180 * tiptimeCessi

    time_labels  = string.(tiptimeticks)
    ax3_top.xticks = (eps_ticks, time_labels)

    xlims!(ax3,     ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])
    xlims!(ax3_top, ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])

    # -------------------------------------------------------------------------
    # Layout tweaks (exactly as original)
    # -------------------------------------------------------------------------
    colgap!(f.layout, 0.0)
    rowgap!(f.layout, 0.0)
    colsize!(f.layout, 2, Relative(0.6))
    rowsize!(f.layout, 2, Relative(0.6))

    f
end

save("./Manuscript/Figures/Fig4.pdf", f, pt_per_unit=1)
#endregion

#region Fig5
@load "./Data/brCessi.jld2"
@load "./Data/brGIS.jld2"

palette = cgrad([safecolor, unsafecolor], 2, categorical = true)

begin
    size_inches = (5.3, 3.7)
    size_pt     = 72 .* size_inches

    f  = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Time series (top-left)
    # -------------------------------------------------------------------------
    ga = f[1, 1] = GridLayout()
    Label(ga[1, 1, TopLeft()], L"(\text{a})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    axm1 = Axis(ga[1, 1:2])

    @load "./Data/TrajGISCessieps1.0.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(axm1, t ./ ϵ .* 180, V; color = nocolor)

    @load "./Data/TrajGISCessieps1.6.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(axm1, t ./ ϵ .* 180, V; color = safecolor)

    @load "./Data/TrajGISCessieps2.4.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(axm1, t ./ ϵ .* 180, V; color = unsafecolor)

    hidexdecorations!(axm1)
    hideydecorations!(axm1, ticks = false, ticklabels = false, label = false)
    ylims!(axm1, 0.0, 1.1)
    xlims!(axm1, 0, 1.0e6)

    axm1.ylabel = L"V"

    ax0 = Axis(ga[2, 1:2])

    @load "./Data/TrajGISCessieps1.0.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax0, t ./ ϵ .* 180, Q; color = nocolor)

    @load "./Data/TrajGISCessieps1.6.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax0, t ./ ϵ .* 180, Q; color = safecolor)

    @load "./Data/TrajGISCessieps2.4.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax0, t ./ ϵ .* 180, Q; color = unsafecolor)

    hideydecorations!(ax0, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax0, ticks = false, ticklabels = false, label = false)

    ylims!(ax0, 0.7, 5.0)
    xlims!(ax0, 0, 1.0e6)

    ax0.ylabel = L"Q"
    ax0.xlabel = L"\text{Model years}"

    rowgap!(ga, 0)
    colgap!(ga, 0)

    # -------------------------------------------------------------------------
    # (b) Bifurcation-style panels (bottom-left)
    # -------------------------------------------------------------------------
    gb = f[2, 1] = GridLayout()
    Label(gb[1, 1, TopLeft()], L"(\text{b})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax1 = Axis(gb[1, 1]; xticks = [1.2, 1.4])

    iend = brGIS.specialpoint[3].idx
    sn1  = brGIS.specialpoint[1].idx
    sn2  = brGIS.specialpoint[2].idx

    lines!(ax1, brGIS.branch.param[1:sn1],   brGIS.branch.V[1:sn1];   color = :black)
    lines!(ax1, brGIS.branch.param[sn1:sn2], brGIS.branch.V[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brGIS.branch.param[sn2:end], brGIS.branch.V[sn2:end]; color = :black)

    @load "./Data/TrajGISCessieps1.0.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :grey)

    xlims!(ax1, 1.1, 2.3)

    ax1.ylabel = L"V"
    ax1.xlabel = L"\delta T"

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax2 = Axis(gb[1, 2], yaxisposition = :right, xticks = [-0.5, -0.25, 0.0])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax2, -1 .* brCessi.branch.param[1:sn1] .+ parametersGISCessi.nupl.F0,   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax2, -1 .* brCessi.branch.param[sn1:sn2] .+ parametersGISCessi.nupl.F0, brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax2, -1 .* brCessi.branch.param[sn2:end] .+ parametersGISCessi.nupl.F0, brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/TrajGISCessieps1.6.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax2, -1 .* γ .* ϵ .* flux, Q; color = safecolor)

    @load "./Data/TrajGISCessieps2.4.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax2, -1 .* γ .* ϵ .* flux, Q; color = unsafecolor)

    @load "./Data/TrajGISCessieps1.0.jld2" r parametersGISCessi ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax2, -1 .* γ .* ϵ .* flux, Q; color = nocolor)

    ax2.xlabel = L"\gamma\dot{V}"
    ax2.ylabel = L"T"

    xlims!(ax2, -0.5, 0.0)

    hideydecorations!(ax2, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax2, ticks = false, ticklabels = false, label = false)

    colgap!(gb, 0)

    # -------------------------------------------------------------------------
    # (c) Parameter plane (right column)
    # -------------------------------------------------------------------------
    gc = f[1:2, 2] = GridLayout()
    Label(gc[1, 1, TopLeft()], L"(\text{c})";
        fontsize = 9,
        font     = :bold,
        padding  = (0, 5, 5, 0),
        halign   = :right
    )

    ax3 = Axis(gc[1, 1])

    @load "./Data/ResultGISCessi.jld2" resultGISCessi ϵrangeGISCessi γrangeGISCessi critGISCessi parametersGISCessi statemax tiptimeGIS critcorr tmax

    contourf!(ax3, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; colormap = palette)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critcorr; levels = [1], color = :black, linestyle = :dash, linewidth = 1)

    band_y = parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = :white)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = nocolor)

    hlines!(ax3, 3.73; color = :grey, linestyle = :dot)
    vlines!(ax3, 0.38; color = :grey, linestyle = :dot)

    ylims!(ax3, 0, 10)

    scatter!(ax3, 2.2, 3.8; marker = :star5, strokewidth = 1, markersize = 7, color = unsafecolor)
    scatter!(ax3, 1.6, 3.8; marker = :star5, strokewidth = 1, markersize = 7, color = safecolor)
    scatter!(ax3, 1.0, 3.8; marker = :star5, strokewidth = 1, markersize = 7, color = nocolor)

    ax3.xlabel = L"Timescale separation $\epsilon$"
    ax3.ylabel = L"Coupling strength $\gamma$"

    # --- Top axis with GIS tipping duration labels ---------------------------
    ax3_top = Axis(gc[1, 1];
        xaxisposition = :top,
        yaxisposition = :right,
        xlabel        = L"\text{GIS tipping duration (years)}",
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [8000, 4000, 2000, 1000]
    eps_ticks    = 1 ./ tiptimeticks * 180 / 470 * tiptimeGIS

    time_labels  = string.(tiptimeticks)
    ax3_top.xticks = (eps_ticks, time_labels)

    xlims!(ax3,     ϵrangeGISCessi[1], ϵrangeGISCessi[end])
    xlims!(ax3_top, ϵrangeGISCessi[1], ϵrangeGISCessi[end])

    # -------------------------------------------------------------------------
    # Layout tweaks (exactly as original)
    # -------------------------------------------------------------------------
    colgap!(f.layout, 0.0)
    rowgap!(f.layout, 0.0)
    colsize!(f.layout, 2, Relative(0.6))
    rowsize!(f.layout, 2, Relative(0.6))

    f
end


save("./Manuscript/Figures/Fig5.pdf", f, pt_per_unit=1)
#endregion

#region FigS1
@load "./Data/brCessi.jld2"
@load "./Data/brVeg.jld2"

palette = cgrad([safecolor, unsafecolor], 2, categorical=true)

begin
    size_inches = (5.3, 5.0)
    size_pt     = 72 .* size_inches

    f = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Row 1
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 1]; xticks = [1.2, 1.4])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax1, brCessi.branch.param[1:sn1],   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax1, brCessi.branch.param[sn1:sn2], brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brCessi.branch.param[sn2:end], brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/ResultCessiVeg.jld2" resultCessiVeg ϵrangeCessiVeg γrangeCessiVeg critCessiVeg parametersCessiVeg tiptimeCessi statemax r ϵ γ state0 tspan prob t x y P T Q
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q; color = leadcolor)
    println(r)
    xlims!(ax1, 1.1, 1.6)

    ax1.ylabel = L"Q"

    text!(ax1, 1.55, 4.6; text=L"r = 10^{-4}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax3 = Axis(f[2, 1])

    contourf!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap = palette, levels = 2)
    contour!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels = [-1], color = :black, linewidth = 1)

    γcrit = findfirst(x -> x > -1, critCessiVeg)

    band_y = parametersCessiVeg.nupl.distance / CascAnalytics.c(statemax, parametersCessiVeg.vec, 0)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = :white)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = nocolor)

    ylims!(ax3, 0.0, 0.88)

    ax3.ylabel = L"Coupling strength $\gamma$"

    text!(ax3, 0.017, 0.83; text=L"r = 10^{-4}", align=[:right, :top], fontsize=7)

    # Top axis for AMOC tipping duration (years)
    ax3_top = Axis(f[2, 1];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [200, 100, 50]
    eps_ticks    = 1 ./ tiptimeticks * 1 / 180 * tiptimeCessi
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])
    xlims!(ax3_top, ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])

    # -------------------------------------------------------------------------
    # (b) Row 2
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 2]; xticks = [1.2, 1.4])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax1, brCessi.branch.param[1:sn1],   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax1, brCessi.branch.param[sn1:sn2], brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brCessi.branch.param[sn2:end], brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/ResultCessiVegr3.jld2" resultCessiVeg ϵrangeCessiVeg γrangeCessiVeg critCessiVeg parametersCessiVeg tiptimeCessi statemax r ϵ γ state0 tspan prob t x y P T Q
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q; color = leadcolor)
    println(r)
    xlims!(ax1, 1.1, 1.6)

    text!(ax1, 1.55, 4.6; text=L"r = 10^{-3}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    hideydecorations!(ax1)

    ax3 = Axis(f[2, 2])

    contourf!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap = palette, levels = 2)
    contour!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels = [-1], color = :black, linewidth = 1)

    γcrit = findfirst(x -> x > -1, critCessiVeg)

    band_y = parametersCessiVeg.nupl.distance / CascAnalytics.c(statemax, parametersCessiVeg.vec, 0)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = :white)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = nocolor)

    ylims!(ax3, 0.0, 0.88)

    ax3.ylabel = L"Coupling strength $\gamma$"

    text!(ax3, 0.017, 0.83; text=L"r = 10^{-3}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax3)

    ax3_top = Axis(f[2, 2];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [200, 100, 50]
    eps_ticks    = 1 ./ tiptimeticks * 1 / 180 * tiptimeCessi
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])
    xlims!(ax3_top, ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])

    # -------------------------------------------------------------------------
    # (c) Row 3
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 3]; xticks = [1.2, 1.4])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax1, brCessi.branch.param[1:sn1],   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax1, brCessi.branch.param[sn1:sn2], brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brCessi.branch.param[sn2:end], brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/ResultCessiVegr2.jld2" resultCessiVeg ϵrangeCessiVeg γrangeCessiVeg critCessiVeg parametersCessiVeg tiptimeCessi statemax r ϵ γ state0 tspan prob t x y P T Q
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q; color = leadcolor)
    println(r)
    xlims!(ax1, 1.1, 1.6)

    text!(ax1, 1.55, 4.6; text=L"r = 10^{-2}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    hideydecorations!(ax1)

    ax3 = Axis(f[2, 3])
    hideydecorations!(ax3)

    contourf!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap = palette, levels = 2)
    contour!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels = [-1], color = :black, linewidth = 1)

    γcrit = findfirst(x -> x > -1, critCessiVeg)

    band_y = parametersCessiVeg.nupl.distance / CascAnalytics.c(statemax, parametersCessiVeg.vec, 0)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = :white)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = nocolor)

    ylims!(ax3, 0.0, 0.88)

    text!(ax3, 0.017, 0.83; text=L"r = 10^{-2}", align=[:right, :top], fontsize=7)


    ax3_top = Axis(f[2, 3];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [200, 100, 50]
    eps_ticks    = 1 ./ tiptimeticks * 1 / 180 * tiptimeCessi
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])
    xlims!(ax3_top, ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])

    # -------------------------------------------------------------------------
    # (c) Row 3
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 4]; xticks = [1.2, 1.4])

    iend = brCessi.specialpoint[3].idx
    sn1  = brCessi.specialpoint[1].idx
    sn2  = brCessi.specialpoint[2].idx

    lines!(ax1, brCessi.branch.param[1:sn1],   brCessi.branch.Q[1:sn1];   color = :black)
    lines!(ax1, brCessi.branch.param[sn1:sn2], brCessi.branch.Q[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brCessi.branch.param[sn2:end], brCessi.branch.Q[sn2:end]; color = :black)

    @load "./Data/ResultCessiVegr1.jld2" resultCessiVeg ϵrangeCessiVeg γrangeCessiVeg critCessiVeg parametersCessiVeg tiptimeCessi statemax r ϵ γ state0 tspan prob t x y P T Q
    lines!(ax1, parametersCessiVeg.nupl.F0 .+ parametersCessiVeg.nupl.r .* t, Q; color = leadcolor)
    println(r)
    xlims!(ax1, 1.1, 1.6)

    text!(ax1, 1.55, 4.6; text=L"r = 10^{-1}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    hideydecorations!(ax1)

    ax3 = Axis(f[2, 4])
    hideydecorations!(ax3)

    contourf!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, resultCessiVeg; colormap = palette, levels = 2)
    contour!(ax3, ϵrangeCessiVeg, γrangeCessiVeg, critCessiVeg; levels = [-1], color = :black, linewidth = 1)

    γcrit = findfirst(x -> x > -1, critCessiVeg)

    band_y = parametersCessiVeg.nupl.distance / CascAnalytics.c(statemax, parametersCessiVeg.vec, 0)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = :white)
    band!(ax3, ϵrangeCessiVeg, zeros(length(ϵrangeCessiVeg)), band_y; color = nocolor)

    ylims!(ax3, 0.0, 0.88)

    text!(ax3, 0.017, 0.83; text=L"r = 10^{-1}", align=[:right, :top], fontsize=7)


    ax3_top = Axis(f[2, 4];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [200, 100, 50]
    eps_ticks    = 1 ./ tiptimeticks * 1 / 180 * tiptimeCessi
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])
    xlims!(ax3_top, ϵrangeCessiVeg[1], ϵrangeCessiVeg[end])

    # -------------------------------------------------------------------------
    # Layout tweaks (keep behavior)
    # -------------------------------------------------------------------------
    colgap!(f.layout, 0.0)
    rowgap!(f.layout, 7.5)

    # colsize!(ga, 1, Relative(0.4))
    # colsize!(gb, 1, Relative(0.4))
    # colsize!(gc, 1, Relative(0.4))

    axdown = Axis(f[2,1:4], xlabelpadding = 17)
    hidedecorations!(axdown, label=false)
    hidespines!(axdown)

    axdown.xlabel = L"\text{Timescale separation } \epsilon"

    axdown_top = Axis(f[2, 1:4];
        xaxisposition = :top,
        yaxisposition = :right,
        xlabel        = L"\text{AMOC tipping duration (years)}",
        xlabelpadding = 17,
    )

    hidedecorations!(axdown_top, grid = false, label = false)
    hidespines!(axdown_top, :l, :r, :b)

    axup = Axis(f[1,1:4], xlabelpadding = 17)
    hidedecorations!(axup, label=false)
    hidespines!(axup)

    axup.xlabel = L"F"

    f
end

save("./Manuscript/Figures/FigS1.pdf", f, pt_per_unit=1)
#endregion

#region FigS2
@load "./Data/brCessi.jld2"
@load "./Data/brGIS.jld2"

palette = cgrad([safecolor, unsafecolor], 2, categorical=true)
CascAnalytics.c(statemax, parametersGISCessi.vec, tmax)

parametersGISCessi.vec
begin
    size_inches = (5.3, 5.0)
    size_pt     = 72 .* size_inches

    f = Figure(size = size_pt, fontsize = 7, padding = 0)

    # -------------------------------------------------------------------------
    # (a) Row 1
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 1]; xticks = [1.3, 1.7, 2.1])

    println(ax1.limits)

    iend = brGIS.specialpoint[3].idx
    sn1  = brGIS.specialpoint[1].idx
    sn2  = brGIS.specialpoint[2].idx

    lines!(ax1, brGIS.branch.param[1:sn1],   brGIS.branch.V[1:sn1];   color = :black)
    lines!(ax1, brGIS.branch.param[sn1:sn2], brGIS.branch.V[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brGIS.branch.param[sn2:end], brGIS.branch.V[sn2:end]; color = :black)

    @load "./Data/ResultGISCessi.jld2" resultGISCessi ϵrangeGISCessi γrangeGISCessi critGISCessi parametersGISCessi statemax tiptimeGIS critcorr tmax r ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :grey)

    ax1.ylabel = L"V"

    xlims!(ax1, 1.1, 2.3)
    ylims!(ax1, 0.0, 1.0)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    text!(ax1, 2.2, 0.95; text=L"r = 10^{-4}", align=[:right, :top], fontsize=7)

    ax3 = Axis(f[2, 1], ylabel = L"Coupling strength $\gamma$")

    contourf!(ax3, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; colormap = palette)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critGISCessi;     levels = [1], color = :black, linewidth = 1)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critcorr;    levels = [1], color = :black, linewidth = 1, linestyle = :dash)

    band_y = parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = :white)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = nocolor)

    ylims!(ax3, 0, 10)

    text!(ax3, 2.7, 9.5; text=L"r = 10^{-4}", align=[:right, :top], fontsize=7)

    ax3_top = Axis(f[2, 1];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [4000, 2000, 1000]
    eps_ticks    = 1 ./ tiptimeticks * 180 / 470 * tiptimeGIS
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeGISCessi[1], ϵrangeGISCessi[end])
    xlims!(ax3_top, ϵrangeGISCessi[1], ϵrangeGISCessi[end])

    # -------------------------------------------------------------------------
    # (b) Row 2
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 2]; xticks = [1.3, 1.7, 2.1])

    iend = brGIS.specialpoint[3].idx
    sn1  = brGIS.specialpoint[1].idx
    sn2  = brGIS.specialpoint[2].idx

    lines!(ax1, brGIS.branch.param[1:sn1],   brGIS.branch.V[1:sn1];   color = :black)
    lines!(ax1, brGIS.branch.param[sn1:sn2], brGIS.branch.V[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brGIS.branch.param[sn2:end], brGIS.branch.V[sn2:end]; color = :black)

    @load "./Data/ResultGISCessir3.jld2" resultGISCessi ϵrangeGISCessi γrangeGISCessi critGISCessi parametersGISCessi statemax tiptimeGIS critcorr tmax r ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :grey)

    ax1.ylabel = L"V"

    xlims!(ax1, 1.1, 2.3)
    ylims!(ax1, 0.0, 1.0)

    hideydecorations!(ax1)

    text!(ax1, 2.2, 0.95; text=L"r = 10^{-3}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax3 = Axis(f[2, 2])

    contourf!(ax3, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; colormap = palette)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critGISCessi;     levels = [1], color = :black, linewidth = 1)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critcorr;    levels = [1], color = :black, linewidth = 1, linestyle = :dash)

    band_y = parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = :white)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = nocolor)

    ylims!(ax3, 0, 10)

    hideydecorations!(ax3, grid = false)

    text!(ax3, 2.7, 9.5; text=L"r = 10^{-3}", align=[:right, :top], fontsize=7)

    ax3_top = Axis(f[2, 2];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [4000, 2000, 1000]
    eps_ticks    = 1 ./ tiptimeticks * 180 / 470 * tiptimeGIS
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeGISCessi[1], ϵrangeGISCessi[end])
    xlims!(ax3_top, ϵrangeGISCessi[1], ϵrangeGISCessi[end])

    # -------------------------------------------------------------------------
    # (c) Row 3
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 3]; xticks = [1.3, 1.7, 2.1])

    iend = brGIS.specialpoint[3].idx
    sn1  = brGIS.specialpoint[1].idx
    sn2  = brGIS.specialpoint[2].idx

    lines!(ax1, brGIS.branch.param[1:sn1],   brGIS.branch.V[1:sn1];   color = :black)
    lines!(ax1, brGIS.branch.param[sn1:sn2], brGIS.branch.V[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brGIS.branch.param[sn2:end], brGIS.branch.V[sn2:end]; color = :black)

    @load "./Data/ResultGISCessir2.jld2" resultGISCessi ϵrangeGISCessi γrangeGISCessi critGISCessi parametersGISCessi statemax tiptimeGIS critcorr tmax r ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :grey)

    ax1.ylabel = L"V"

    xlims!(ax1, 1.1, 2.3)
    ylims!(ax1, 0.0, 1.0)

    hideydecorations!(ax1)

    text!(ax1, 2.2, 0.95; text=L"r = 10^{-2}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax3 = Axis(f[2, 3])

    contourf!(ax3, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; colormap = palette)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critGISCessi;     levels = [1], color = :black, linewidth = 1)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critcorr;    levels = [1], color = :black, linewidth = 1, linestyle = :dash)

    band_y = parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = :white)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = nocolor)

    ylims!(ax3, 0, 10)

    hideydecorations!(ax3, grid = false)

    text!(ax3, 2.7, 9.5; text=L"r = 10^{-2}", align=[:right, :top], fontsize=7)

    ax3_top = Axis(f[2, 3];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [2000, 1000]
    eps_ticks    = 1 ./ tiptimeticks * 180 / 470 * tiptimeGIS
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeGISCessi[1], ϵrangeGISCessi[end])
    xlims!(ax3_top, ϵrangeGISCessi[1], ϵrangeGISCessi[end])

    # -------------------------------------------------------------------------
    # (c) Row 3
    # -------------------------------------------------------------------------

    ax1 = Axis(f[1, 4]; xticks = [1.3, 1.7, 2.1])

    iend = brGIS.specialpoint[3].idx
    sn1  = brGIS.specialpoint[1].idx
    sn2  = brGIS.specialpoint[2].idx

    lines!(ax1, brGIS.branch.param[1:sn1],   brGIS.branch.V[1:sn1];   color = :black)
    lines!(ax1, brGIS.branch.param[sn1:sn2], brGIS.branch.V[sn1:sn2]; color = :black, linestyle = :dash)
    lines!(ax1, brGIS.branch.param[sn2:end], brGIS.branch.V[sn2:end]; color = :black)

    @load "./Data/ResultGISCessir1.jld2" resultGISCessi ϵrangeGISCessi γrangeGISCessi critGISCessi parametersGISCessi statemax tiptimeGIS critcorr tmax r ϵ γ state0 tspan prob t V x y Q flux
    lines!(ax1, 1.1 .+ parametersGISCessi.nupl.r .* t, V; color = :grey)

    ax1.ylabel = L"V"

    xlims!(ax1, 1.1, 2.3)
    ylims!(ax1, 0.0, 1.0)

    hideydecorations!(ax1)

    text!(ax1, 2.2, 0.95; text=L"r = 10^{-1}", align=[:right, :top], fontsize=7)

    hideydecorations!(ax1, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax1, ticks = false, ticklabels = false, label = false)

    ax3 = Axis(f[2, 4])

    contourf!(ax3, ϵrangeGISCessi, γrangeGISCessi, resultGISCessi; colormap = palette)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critGISCessi;     levels = [1], color = :black, linewidth = 1)
    contour!(ax3, ϵrangeGISCessi, γrangeGISCessi, critcorr;    levels = [1], color = :black, linewidth = 1, linestyle = :dash)

    band_y = parametersGISCessi.nupl.distance ./ (CascAnalytics.c(statemax, parametersGISCessi.vec, tmax) .* ϵrangeGISCessi)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = :white)
    band!(ax3, ϵrangeGISCessi, zeros(length(ϵrangeGISCessi)), band_y; color = nocolor)

    ylims!(ax3, 0, 10)

    hideydecorations!(ax3, grid = false)

    text!(ax3, 2.7, 9.5; text=L"r = 10^{-1}", align=[:right, :top], fontsize=7)

    ax3_top = Axis(f[2, 4];
        xaxisposition = :top,
        yaxisposition = :right,
    )

    hideydecorations!(ax3_top, grid = false)
    hidespines!(ax3_top, :l, :r, :b)

    tiptimeticks = [2000, 1000]
    eps_ticks    = 1 ./ tiptimeticks * 180 / 470 * tiptimeGIS
    ax3_top.xticks = (eps_ticks, string.(tiptimeticks))

    xlims!(ax3,     ϵrangeGISCessi[1], ϵrangeGISCessi[end])
    xlims!(ax3_top, ϵrangeGISCessi[1], ϵrangeGISCessi[end])

    # -------------------------------------------------------------------------
    # Layout tweaks 
    # -------------------------------------------------------------------------
    colgap!(f.layout, 0.0)
    rowgap!(f.layout, 7.5)

    # colsize!(ga, 1, Relative(0.4))
    # colsize!(gb, 1, Relative(0.4))
    # colsize!(gc, 1, Relative(0.4))

    axdown = Axis(f[2,1:4], xlabelpadding = 17)
    hidedecorations!(axdown, label=false)
    hidespines!(axdown)

    axdown.xlabel = L"\text{Timescale separation } \epsilon"

    axdown_top = Axis(f[2, 1:4];
        xaxisposition = :top,
        yaxisposition = :right,
        xlabel        = L"\text{GIS tipping duration (years)}",
        xlabelpadding = 17,
    )

    hidedecorations!(axdown_top, grid = false, label = false)
    hidespines!(axdown_top, :l, :r, :b)

    axup = Axis(f[1,1:4], xlabelpadding = 17)
    hidedecorations!(axup, label=false)
    hidespines!(axup)

    axup.xlabel = L"\delta T"

    f
end

save("./Manuscript/Figures/FigS2.pdf", f, pt_per_unit=1)
#endregion