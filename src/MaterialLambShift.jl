module MaterialLambShift

#using Bessels #Does not support complex arguments
using SpecialFunctions
using Plots
using LinearAlgebra
using StaticArrays
using QuadGK
using Debugger
using HCubature

# For proper execution by pyjulia exported names must be ASCII #
#                           :c                                 #

# export dlayerresponse
export dlayerradialresponse
export dlayerradialdynresponse
export hankelplot
export hankelfcnvalues
export dlayerresponsecoef
export coeftoresponse
export staticresponse
export dynamicresponse
export diffresponse
export asymptdynamicresponse
export kcheck
export qcheck
export localresponseωplot
export localresponseiωplot
export disccheck
export localresponse
export asymptcheck
export getν
export ΔEdrudeint1
export ΔEdrudeint1old
export ΔEdrudeχ̃
export ΔEdrudeχ̃plot
export ΔEdrudeisoold
export ΔEdrudecoefold
export ΔEdrudeintold
export drudeτplot
export druder0plot
export mybesselj

## Variable setup
const rw::Float64 = 0.1 # Wire width, ls
const iϵ::ComplexF64 = 1e-11im # Causality fixer

# Returns tuple of Hankel function coordinates (αE, βE, γE, αB, βB, γB)
# consistent with a fluctuating dipole at distance rd from the axis of
# symmetry.
function dlayerresponsecoef(ω, m, k, d; rd=1) 
    ν = √Complex(ω^2 - k^2)
    # Save Hankel fcn values
    H1w, H2w, H′1w, H′2w, H3w, H′3w = hankelfcnvalues(rw, m, ν)
    H1d, H2d, H′1d, H′2d, H3d, H′3d = hankelfcnvalues(rd, m, ν)
    #println((H1w, H2w, H′1w, H′2w, H3w, H′3w))

    # Construct boundary condition matrix
    # (same form as in overleaf document)
    BCmat = [-rd^2*ν*H′1d; H1w; -H1d;     0;           0;              0;;
             -rd^2*ν*H′2d; H2w; -H2d;     0;           0;              0;;
              rd^2*ν*H′3d;   0;  H3d;     0;           0;              0;;
                      0;   0;    0; -H′1d;        -rd*H1d;           H′1w;;
                      0;   0;    0; -H′2d;        -rd*H2d;           H′2w;;
                      0;   0;    0;  H′3d;         rd*H3d;              0;;]
    BCvec = [k*m*d[2] - rd*ν^2*d[3]; 0; 0; 0; im*ω*d[2]; 0]

    return BCmat\BCvec
end # function dlayerresponsecoef
# Returns field values at r due to coef α, by the LSE solution
function coeftoresponse(r, rd, ω, m, k, α)

    ν = √Complex(ω^2 - k^2)
    αE, βE, γE, αB, βB, γB = α
    # Retrieve Hankel fcn values
    H1, H2, H′1, H′2, H3, H′3 = hankelfcnvalues(r, m, ν)
    #println(α)
    #println(k/ν)

    if r < rd
        Er = im*k/ν*(αE*H′1 + βE*H′2) - ω*m/(r*ν^2)*(αB*H1 + βB*H2)
        Eθ = -k*m/(r*ν^2)*(αE*H1 + βE*H2) - im*ω/ν*(αB*H′1 + βB*H′2)
        Ez = αE*H1 +  βE*H2
        Br = m*ω/(r*ν^2)*(αE*H1 + βE*H2) + im*k/ν*(αB*H′1 + βB*H′2)
        Bθ = im*ω/ν*(αE*H′1 + βE*H′2) - k*m/(r*ν^2)*(αB*H1 + βB*H2)
        Bz = αB*H1 + βB*H2
    else
        Er = im*k/ν*γE*H′3 - ω*m/(r*ν^2)*γB*H3
        Eθ = -k*m/(r*ν^2)*γE*H3 - im*ω/ν*γB*H′3
        Ez = γE*H3
        Br = m*ω/(r*ν^2)*γE*H3 + im*k/ν*γB*H′3
        Bθ = im*ω/ν*γE*H′3 - k*m/(r*ν^2)*γB*H3
        Bz = γB*H3
    end # if
    return Er, Eθ, Ez, Br, Bθ, Bz
end # function coeftoresponse

# Returns values of the Hankel functions at point νr as well as their first
# derivatives. If ν is pure imaginary modified Bessel values are returned
# instead.
function hankelfcnvalues(r, m, ν)
    if angle(ν) == 0
        ν = real(ν)
        H1 = hankelh1(m, ν*r)
        H2 = hankelh2(m, ν*r)
        H′1 = 1/2*(hankelh1(m-1, ν*r) - hankelh1(m+1, ν*r))
        H′2 = 1/2*(hankelh2(m-1, ν*r) - hankelh2(m+1, ν*r))
        return (H1, H2, H′1, H′2, H1, H′1) # First species repeated outside of dipole layer
    elseif angle(ν) == pi/2
        ν = imag(ν)
        I = besseli(m, ν*r)
        K = besselk(m, ν*r)
        # Derivatives are multiplied with -im to compensate for
        # d/dr(H)=νH′ in all derivations
        I′ = -im/2*(besseli(m-1, ν*r) + besseli(m+1, ν*r))
        K′ = im/2*(besselk(m-1, ν*r) + besselk(m+1, ν*r))
        return (I, K, I′, K′, K, K′) # Second species repeated
    else
        throw(error())
    end
end # function hankelfcnvalues

# Evaluates electric field at a point (r, θ, z) due to a dipole layer with
# (ω, m, k, rd) NOT UPDATED
function dlayerresponse(r, θ, z, ω, m, k)
    ν = √Complex(ω^2 - k^2)
    # Retrieve coefficients
    αE, βE, γE, αB, βB, γB = dlayerresponsecoef(ω, m, k, d, rd)

    H1, H2, H′1, H′2, H3, H′3 = hankelfcnvalues(r, m, ν)

    if r < rd
        return exp(im*(m*θ + k*z))*(
                                    [im*k/ν*(αE*H′1 + βE*H′2) - ω*m/(r*ν^2)*(αB*H1 + βB*H2);
                                     -k*m/(r*ν^2)*(αE*H1 + βE*H2) - im*ω/ν*(αB*H′1 + βB*H′2);
                                     αE*H1 +  βE*H2;
                                     ω*m/(r*ν^2)*(αE*H1 + βE*H2) + im*k/ν*(αB*H′1 + βB*H′2);
                                     im*ω/ν*(αE*H′1 + βE*H′2) - k*m/(r*ν^2)*(αB*H1 + βB*H2);
                                     αB*H1 + βB*H2]
                                   );
    else return exp(im*(m*θ + k*z))*(
                                    [im*k/ν*γE*H′3 - ω*m/(r*ν^2)*γB*H3;
                                     -k*m/(r*ν^2)*γE*H3 - im*ω/ν*γB*H′3;
                                     γE*H2;
                                     ω*m/(r*ν^2)*γE*H2 + im*k/ν*γB*H′2;
                                     im*ω/ν*γE*H′3 - k*m/(r*ν^2)*γB*H3;
                                     γB*H3]
                                   );
    end # if
end # function dlayerresponse

# Plots dynamic response as a function of r
function dlayerradialdynresponse(ω, m, k; rd=1, maxr=2, stepr=0.001, minr=0.1)
    rr = minr:stepr:maxr # Range of r
    xrr, xrθ, xrz, xθr, xθθ, xθz, xzr, xzθ, xzz = [[] for _ in 1:9] # Dynamic solution
    drr, drθ, drz, dθr, dθθ, dθz, dzr, dzθ, dzz = [[] for _ in 1:9] # Difference solution
    arr, arθ, arz, aθr, aθθ, aθz, azr, azθ, azz = [[] for _ in 1:9] # Asympt. solution
    brr, brθ, brz, bθr, bθθ, bθz, bzr, bzθ, bzz = [[] for _ in 1:9] # Asympt. diff. solution

    for r in rr
        # Dynamic response
        χ = dynamicresponse(r, rd, ω + iϵ, m, k)
        push!(xrr, χ[1, 1])
        push!(xrθ, χ[1, 2])
        push!(xrz, χ[1, 3])
        push!(xθr, χ[2, 1])
        push!(xθθ, χ[2, 2])
        push!(xθz, χ[2, 3])
        push!(xzr, χ[3, 1])
        push!(xzθ, χ[3, 2])
        push!(xzz, χ[3, 3])
    end # for

    for r in rr
        # Diff response
        χ = diffresponse(r, rd, ω + iϵ, m, k)
        push!(drr, χ[1, 1])
        push!(drθ, χ[1, 2])
        push!(drz, χ[1, 3])
        push!(dθr, χ[2, 1])
        push!(dθθ, χ[2, 2])
        push!(dθz, χ[2, 3])
        push!(dzr, χ[3, 1])
        push!(dzθ, χ[3, 2])
        push!(dzz, χ[3, 3])
    end # for

    for r in rr
        # Asymptotic response
        χ = asymptdynamicresponse(r, rd, ω + iϵ, m, k)
        push!(arr, χ[1, 1])
        push!(arθ, χ[1, 2])
        push!(arz, χ[1, 3])
        push!(aθr, χ[2, 1])
        push!(aθθ, χ[2, 2])
        push!(aθz, χ[2, 3])
        push!(azr, χ[3, 1])
        push!(azθ, χ[3, 2])
        push!(azz, χ[3, 3])
    end # for

    for r in rr
        # Asymptotic diff. response
        χ = asymptdiffresponse(r, rd, ω + iϵ, m, k)
        push!(brr, χ[1, 1])
        push!(brθ, χ[1, 2])
        push!(brz, χ[1, 3])
        push!(bθr, χ[2, 1])
        push!(bθθ, χ[2, 2])
        push!(bθz, χ[2, 3])
        push!(bzr, χ[3, 1])
        push!(bzθ, χ[3, 2])
        push!(bzz, χ[3, 3])
    end # for

    # Mask low values, for now
    ymin = 1e-8
    xrr[abs.(xrr) .< ymin] .= ymin
    xrθ[abs.(xrθ) .< ymin] .= ymin
    xrz[abs.(xrz) .< ymin] .= ymin
    xθr[abs.(xθr) .< ymin] .= ymin
    xθθ[abs.(xθθ) .< ymin] .= ymin
    xθz[abs.(xθz) .< ymin] .= ymin
    xzr[abs.(xzr) .< ymin] .= ymin
    xzθ[abs.(xzθ) .< ymin] .= ymin
    xzz[abs.(xzz) .< ymin] .= ymin
    drr[abs.(drr) .< ymin] .= ymin
    drθ[abs.(drθ) .< ymin] .= ymin
    drz[abs.(drz) .< ymin] .= ymin
    dθr[abs.(dθr) .< ymin] .= ymin
    dθθ[abs.(dθθ) .< ymin] .= ymin
    dθz[abs.(dθz) .< ymin] .= ymin
    dzr[abs.(dzr) .< ymin] .= ymin
    dzθ[abs.(dzθ) .< ymin] .= ymin
    dzz[abs.(dzz) .< ymin] .= ymin
    arr[abs.(arr) .< ymin] .= ymin
    arθ[abs.(arθ) .< ymin] .= ymin
    arz[abs.(arz) .< ymin] .= ymin
    aθr[abs.(aθr) .< ymin] .= ymin
    aθθ[abs.(aθθ) .< ymin] .= ymin
    aθz[abs.(aθz) .< ymin] .= ymin
    azr[abs.(azr) .< ymin] .= ymin
    azθ[abs.(azθ) .< ymin] .= ymin
    azz[abs.(azz) .< ymin] .= ymin
    brr[abs.(brr) .< ymin] .= ymin
    brθ[abs.(brθ) .< ymin] .= ymin
    brz[abs.(brz) .< ymin] .= ymin
    bθr[abs.(bθr) .< ymin] .= ymin
    bθθ[abs.(bθθ) .< ymin] .= ymin
    bθz[abs.(bθz) .< ymin] .= ymin
    bzr[abs.(bzr) .< ymin] .= ymin
    bzθ[abs.(bzθ) .< ymin] .= ymin
    bzz[abs.(bzz) .< ymin] .= ymin

    Errplot = plot()
    #plot!(rr, abs.(xrr), line_z=imag(log.(xrr)), label="χ^{rr}")
    #plot!(rr, abs.(arr), line_z=imag(log.(arr)), label="asympt. χ^{rr}")
    #plot!(rr, abs.(brr), line_z=angle.(brr), label="asympt. diff. χ^{rr}")
    plot!(rr, abs.(drr), line_z=angle.(drr), label="diff. χ^{rr}", linewidth=0.5)
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(drr)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xrr./arr), color=:purple)
    Erθplot = plot()
    #plot!(rr, abs.(xrθ), line_z=imag(log.(xrθ)), label="χ^{rθ}")
    plot!(rr, abs.(drθ), line_z=angle.(drθ), label="diff. χ^{rθ}", linewidth=0.5)
    #plot!(rr, abs.(arθ), line_z=imag(log.(arθ)), label="asympt. χ^{rθ}")
    #plot!(rr, abs.(brθ), line_z=angle.(brθ), label="asympt. diff. χ^{rθ}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(drθ)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xrθ./arθ), color=:purple)
    Erzplot = plot()
    #plot!(rr, abs.(xrz), line_z=imag(log.(xrz)), label="χ^{rz}")
    plot!(rr, abs.(drz), line_z=angle.(drz), label="diff. χ^{rz}", linewidth=0.5)
    #plot!(rr, abs.(arz), line_z=imag(log.(arz)), label="asympt. χ^{rz}")
    #plot!(rr, abs.(brz), line_z=angle.(brz), label="asympt. diff. χ^{rz}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(drz)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xrz./arz), color=:purple)
    Eθrplot = plot()
    #plot!(rr, abs.(xθr), line_z=imag(log.(xθr)), label="χ^{θr}")
    plot!(rr, abs.(dθr), line_z=angle.(dθr), label="diff. χ^{θr}", linewidth=0.5)
    #plot!(rr, abs.(aθr), line_z=imag(log.(aθr)), label="asympt. χ^{θr}")
    #plot!(rr, abs.(bθr), line_z=angle.(bθr), label="asympt. diff. χ^{θr}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dθr)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xθr./aθr), color=:purple)
    Eθθplot = plot()
    #plot!(rr, abs.(xθθ), line_z=imag(log.(xθθ)), label="χ^{θθ}")
    plot!(rr, abs.(dθθ), line_z=angle.(dθθ), label="diff. χ^{θθ}", linewidth=0.5)
    #plot!(rr, abs.(aθθ), line_z=imag(log.(aθθ)), label="asympt. χ^{θθ}")
    #plot!(rr, abs.(bθθ), line_z=angle.(bθθ), label="asympt. diff. χ^{θθ}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dθθ)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xθθ./aθθ), color=:purple)
    Eθzplot = plot()
    #plot!(rr, abs.(xθz), line_z=imag(log.(xθz)), label="χ^{θz}")
    plot!(rr, abs.(dθz), line_z=angle.(dθz), label="diff. χ^{θz}", linewidth=0.5)
    #plot!(rr, abs.(aθz), line_z=imag(log.(aθz)), label="asympt. χ^{θz}")
    #plot!(rr, abs.(bθz), line_z=angle.(bθz), label="asympt. diff. χ^{θz}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dθz)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xθz./aθz), color=:purple)
    Ezrplot = plot()
    #plot!(rr, abs.(xzr), line_z=imag(log.(xzr)), label="χ^{zr}")
    plot!(rr, abs.(dzr), line_z=angle.(dzr), label="diff. χ^{zr}", linewidth=0.5)
    #plot!(rr, abs.(azr), line_z=imag(log.(azr)), label="asympt. χ^{zr}")
    #plot!(rr, abs.(bzr), line_z=angle.(bzr), label="asympt. diff. χ^{zr}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dzr)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xzr./azr), color=:purple)
    Ezθplot = plot()
    #plot!(rr, abs.(xzθ), line_z=imag(log.(xzθ)), label="χ^{zθ}")
    plot!(rr, abs.(dzθ), line_z=angle.(dzθ), label="diff. χ^{zθ}", linewidth=0.5)
    #plot!(rr, abs.(azθ), line_z=imag(log.(azθ)), label="asympt. χ^{zθ}")
    #plot!(rr, abs.(bzθ), line_z=angle.(bzθ), label="asympt. diff. χ^{zθ}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dzθ)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xzθ./azθ), color=:purple)
    Ezzplot = plot()
    #plot!(rr, abs.(xzz), line_z=imag(log.(xzz)), label="χ^{zz}")
    plot!(rr, abs.(dzz), line_z=angle.(dzz), label="diff. χ^{zz}", linewidth=0.5)
    #plot!(rr, abs.(azz), line_z=imag(log.(azz)), label="asympt. χ^{zz}")
    #plot!(rr, abs.(bzz), line_z=angle.(bzz), label="asympt. diff. χ^{zz}")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(dzz)), ("r'", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(xzz./azz), color=:purple)
    plot(Errplot, Erθplot, Erzplot, Eθrplot, Eθθplot, Eθzplot, Ezrplot, Ezθplot, Ezzplot, 
         layout=(3,3), yscale=:log10, plot_title="ω = $ω, m = $m, k = $k", 
         legend_position=:bottomright, legend_font_pointsize=3, size=(800, 500), clims=(0, pi))

end # function dlayerradialdynresponse
# Plots response as a function of r, currently compares dynamic and LSE
# solution
function dlayerradialresponse(ω, m, k, d=[0; 1; 0]; rd=1, maxr=2, stepr=0.001, minr=0.1)
    # Retrieve coefficients
    α = dlayerresponsecoef(ω, m, k, d, rd=rd)

    rr = minr:stepr:maxr # Range of r
    fr = [] # Corresponding responses, E-field
    fθ = [] 
    fz = [] 
    gr = [] # B-field
    gθ = [] 
    gz = [] 
    sr = [] # Static solution
    sθ = []
    sz = []
    xr = [] # Dynamic solution
    xθ = []
    xz = []


    for r in rr
        Er, Eθ, Ez, Br, Bθ, Bz = coeftoresponse(r, rd, ω, m, k, α)./(2*pi)^3 # Factor due to dirac delta decomposition
        # Store in fs, gs
        push!(fr, Er)
        push!(fθ, Eθ)
        push!(fz, Ez)
        push!(gr, Br)
        push!(gθ, Bθ)
        push!(gz, Bz)
        # Static response for comparison
        χ = staticresponse(r, rd, m, k)
        push!(sr, (χ*d)[1])
        push!(sθ, (χ*d)[2])
        push!(sz, (χ*d)[3])
        # Dynamic response for comparison
        χ = dynamicresponse(r, rd, ω + iϵ, m, k)
        push!(xr, (χ*d)[1])
        push!(xθ, (χ*d)[2])
        push!(xz, (χ*d)[3])

    end # for
    # Mask low values, for now
    ymin = 1e-5
    fr[abs.(fr) .< ymin] .= ymin
    gr[abs.(gr) .< ymin] .= ymin
    fθ[abs.(fθ) .< ymin] .= ymin
    gθ[abs.(gθ) .< ymin] .= ymin
    fz[abs.(fz) .< ymin] .= ymin
    gz[abs.(gz) .< ymin] .= ymin
    sr[abs.(sr) .< ymin] .= ymin
    sθ[abs.(sθ) .< ymin] .= ymin
    sz[abs.(sz) .< ymin] .= ymin
    xr[abs.(xr) .< ymin] .= ymin
    xθ[abs.(xθ) .< ymin] .= ymin
    xz[abs.(xz) .< ymin] .= ymin
    #println("sr/xr = $((sr./xr)[20])")
    Erplot = plot()
    plot!(rr, abs.(fr), line_z=angle.(fr), label="E_r")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(fr)), ("r_d", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(sr), #line_z=angle.(sr)), 
          #linestyle=:dash, label="Static E_r")
#    plot!(rr, abs.(xr), line_z=angle.(xr)), 
#          linestyle=:dash, label="Dynamic E_r")
    #plot!(rr, abs.(sr./xr), color=:purple)
    Brplot = plot()
    plot!(rr, abs.(gr), line_z=angle.(gr), label="B_r")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(gr)), ("r_d", 6, :bottom, :left, :red)))
    Eθplot = plot()
    plot!(rr, abs.(fθ), line_z=angle.(fθ), label="E_θ")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(fθ)), ("r_d", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(sθ),line_z=angle.(sθ)), 
          #linestyle=:dash, label="Static E_θ")
    #plot!(rr, abs.(xθ), line_z=angle.(xθ)), 
          #linestyle=:dash, label="Dynamic E_θ")
    #plot!(rr, abs.(sθ./xθ), color=:purple)
    Bθplot = plot(rr, abs.(gθ), line_z=angle.(gθ), label="B_θ")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(gθ)), ("r_d", 6, :bottom, :left, :red)))
    Ezplot = plot()
    plot!(rr, abs.(fz), line_z=angle.(fz), label="E_z")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(fz)), ("r_d", 6, :bottom, :left, :red)))
    #plot!(rr, abs.(sz),# line_z=angle.(sz)), 
          #linestyle=:dash, label="Static E_z")
    #plot!(rr, abs.(xz), line_z=angle.(xz)), 
          #linestyle=:dash, label="Dynamic E_z")
    #plot!(rr, abs.(sz./xz), color=:purple)
    Bzplot = plot(rr, abs.(gz), line_z=angle.(gz), label="B_z")
    plot!([rd], seriestype=:vline, label="", color=:red, annotations=(rd*1.1, minimum(abs.(gz)), ("r_d", 6, :bottom, :left, :red)))
    plot(Erplot, Brplot, Eθplot, Bθplot, Ezplot, Bzplot, layout=(3,2), yscale=:log10, plot_title="ω = $ω, m = $m, k = $k, dz = $(d[3]), dθ=$(d[2]), dr = $(d[1])", legend_position=:topright, legend_font_pointsize=4, clims=(-pi, pi))

end # function dlayerradialresponse

# Returns static response (Er, Eθ, Ez)^t * (dr, dθ, dz) at distance r, r′ for m, k
function staticresponse(r, r′, m, k)
    if r > r′
        return adjoint(staticresponse(r′, r, m, k))
    else
        ka = abs(k)
        #Retrieve  mod. Bessel fcn values
        Iw = besseli(m, ka*rw)
        Kw = besselk(m, ka*rw)
        Ir = besseli(m, ka*r)
        Kr = besselk(m, ka*r)
        Kr′ = besselk(m, ka*r′)
        I′r = 1/2*(besseli(m-1, ka*r) + besseli(m+1, ka*r))
        K′r = -1/2*(besselk(m-1, ka*r) + besselk(m+1, ka*r))
        K′r′ = -1/2*(besselk(m-1, ka*r′) + besselk(m+1, ka*r′))
    
        Qr = Iw/Kw*Kr - Ir
        Q′r = Iw/Kw*K′r - I′r
    
        # Response if r < r′
        χ = (2*pi)^(-2) * [ka*Q′r; im*m*Qr/r; im*k*Qr] * [ka*K′r′;; -im*m*Kr′/r′;; -im*k*Kr′]
        return χ ./(2*pi)
    end
    
end # function staticresponse

# Returns dynamical response(Er, Eθ, Ez)^t * (dr, dθ, dz) at distance r,
# r′ for ω, m, k
function dynamicresponse(r, r′, ω, m, k; verbatim=false)
    kcfac = 1e2 #factor to determine kc
    ν = getν(ω, k)
    # Fix large k errors by returning 0 instead
    if abs(k) > abs(ω) && abs(ν) > kcfac/abs(r-r′)
        return zeros(Complex, 3, 3)
    end
    #Retrieve Hankel fcn values
    H1w = hankelh1(m, ν*rw) #
    H2w = hankelh2(m, ν*rw) #
    H′1w = 1/2*(hankelh1(m-1, ν*rw) - hankelh1(m+1, ν*rw)) #
    H′2w = 1/2*(hankelh2(m-1, ν*rw) - hankelh2(m+1, ν*rw)) #
    H1 = hankelh1(m, ν*r) #
    H1′ = hankelh1(m, ν*r′)
    H′1 = 1/2*(hankelh1(m-1, ν*r) - hankelh1(m+1, ν*r)) #
    H′1′ = 1/2*(hankelh1(m-1, ν*r′) - hankelh1(m+1, ν*r′))  #
#    if verbatim
#        println("ν = $ν")
#        println("H1w = $H1w")
#        println("H2w = $H2w")
#        println("H′1w = $H′1w")
#        println("H′2w = $H′2w")
#        println("H1 = $H1")
#        println("H2 = $H2")
#        println("H′1 = $H′1")
#        println("H′2 = $H′2")
#        println("H1′ = $H1′")
#        println("H2′ = $H2′")
#        println("H′1′ = $H′1′")
#        println("H′2′ = $H′2′")
#        println("Qr = $Qr")
#        println("Qr′ = $Qr′")
#        println("QBr = $QBr")
#        println("QBr′ = $QBr′")
#        println("Q′r = $Q′r")
#        println("Q′r′ = $Q′r′")
#        println("QB′r = $QB′r")
#        println("QB′r′ = $QB′r′")
#    end # if
    if r > r′
        #Retrieve Hankel fcn values
        H′2′ = 1/2*(hankelh2(m-1, ν*r′) - hankelh2(m+1, ν*r′)) 
        H2′ = hankelh2(m, ν*r′)
        Qr′ = H2′ - H2w/H1w*H1′ 
        QBr′ = H2′ - H′2w/H′1w*H1′ 
        Q′r′ = H′2′ - H2w/H1w*H′1′ 
        QB′r′ = H′2′ - H′2w/H′1w*H′1′ 
        # Part due to G_E
        return im/(32*pi^2)*(
                             [im*k*H′1; -k*m/(r*ν)*H1; ν*H1] *
                             [-im*k*Q′r′;; -k*m/(ν*r′)*Qr′;; ν*Qr′]
                             + ω^2*[-im*m/(r*ν)*H1; H′1; 0] *
                             [im*m/(r′*ν)*QBr′;; QB′r′;; 0])
    else
        H2 = hankelh2(m, ν*r)
        H′2 = 1/2*(hankelh2(m-1, ν*r) - hankelh2(m+1, ν*r))
        Qr = H2 - H2w/H1w*H1
        QBr = H2 - H′2w/H′1w*H1
        Q′r = H′2 - H2w/H1w*H′1
        QB′r = H′2 - H′2w/H′1w*H′1
        return im/(32*pi^2)*(
                             [im*k*Q′r; -k*m/(r*ν)*Qr; ν*Qr] *
                             [-im*k*H′1′;; -k*m/(ν*r′)*H1′;; ν*H1′]
                             + ω^2*[-im*m/(r*ν)*QBr; QB′r; 0] *
                             [im*m/(r′*ν)*H1′;; H′1′;; 0])
    end
end # function dynamicresponse

# Returns difference of dynamic wire and free response (Er, Eθ, Ez)^t * (dr, dθ, dz) at distance r,
# r′ for ω, m, k
# Has to be called with causality fix +iϵ to ω
function diffresponse(r, r′, ω, m, k; verbatim=false, diagonal=false)
    # kcfac = 1e2 #factor to determine kc
    ν = getν(ω, k)
    # Fix large k errors by returning 0 instead
#    if abs(k) > abs(ω) && abs(ν) > kcfac/abs(r-r′)
#        return zeros(Complex, 3, 3)
#    end

    # Careful here, should increase performance but might set some values to 0
    #set_zero_subnormals(true) 
    
    Jw = mybesselj(m, ν*rw)
    Yw = mybessely(m, ν*rw)
    J′w = 1/2*(mybesselj(m-1, ν*rw) - mybesselj(m+1, ν*rw))
    Y′w = 1/2*(mybessely(m-1, ν*rw) - mybessely(m+1, ν*rw))
    Qpref = -2*Jw/(Jw + im*Yw)
    QBpref = -2*J′w/(J′w + im*Y′w)
    
    # I suspect that the nan-check is no longer needed anymore, couldn't find
    # the right parameters to test it though
    if isnan(Qpref)
        Qpref = 0
    end # if
    if isnan(QBpref)
        QBpref = 0
    end # if
    H1 = hankelh1(m, ν*r) #
    H1′ = hankelh1(m, ν*r′)
    H′1 = 1/2*(hankelh1(m-1, ν*r) - hankelh1(m+1, ν*r)) #
    H′1′ = 1/2*(hankelh1(m-1, ν*r′) - hankelh1(m+1, ν*r′))  #
#    if verbatim
#        println("ν = $ν")
#        println("H1w = $H1w")
#        println("H2w = $H2w")
#        println("H′1w = $H′1w")
#        println("H′2w = $H′2w")
#        println("H1 = $H1")
#        println("H2 = $H2")
#        println("H′1 = $H′1")
#        println("H′2 = $H′2")
#        println("H1′ = $H1′")
#        println("H2′ = $H2′")
#        println("H′1′ = $H′1′")
#        println("Jw = $Jw")
#        println("Yw = $Yw")
#        println("J′w = $J′w")
#        println("Y′w = $Y′w")
#        println("H′2′ = $H′2′")
#        println("Qr = $Qr")
#        println("Qr′ = $Qr′")
#        println("QBr = $QBr")
#        println("QBr′ = $QBr′")
#        println("Q′r = $Q′r")
#        println("Q′r′ = $Q′r′")
#        println("QB′r = $QB′r")
#        println("QB′r′ = $QB′r′")
#    end # if
    if !diagonal
        if r > r′
            #Retrieve Hankel fcn values
            Qr′ = Qpref*H1′
            QBr′ = QBpref*H1′ 
            Q′r′ = Qpref*H′1′ 
            QB′r′ = QBpref*H′1′ 
            χ = SMatrix{3, 3}(im/(32*pi^2)*(
                                 [im*k*H′1; -k*m/(r*ν)*H1; ν*H1] *
                                 [-im*k*Q′r′;; -k*m/(ν*r′)*Qr′;; ν*Qr′]
                                 + ω^2*[-im*m/(r*ν)*H1; H′1; 0] *
                                 [im*m/(r′*ν)*QBr′;; QB′r′;; 0]))
        else
            Qr = Qpref*H1
            QBr = QBpref*H1 
            Q′r = Qpref*H′1 
            QB′r = QBpref*H′1 
            χ = SMatrix{3, 3}(im/(32*pi^2)*(
                                 [im*k*Q′r; -k*m/(r*ν)*Qr; ν*Qr] *
                                 [-im*k*H′1′;; -k*m/(ν*r′)*H1′;; ν*H1′]
                                 + ω^2*[-im*m/(r*ν)*QBr; QB′r; 0] *
                                 [im*m/(r′*ν)*H1′;; H′1′;; 0]))
        end
    else # Avoids integrating unnecessary elements for localresponse
        if r > r′
            Qr′ = Qpref*H1′
            QBr′ = QBpref*H1′ 
            Q′r′ = Qpref*H′1′ 
            QB′r′ = QBpref*H′1′ 
            χ = SMatrix{3, 3}(im/(32*pi^2)*diagm([k^2*H′1*Q′r′ + ω^2*m^2/(r*r′*ν^2)*H1*QBr′,
                                        k^2*m^2/(ν^2*r*r′)*H1*Qr′ + ω^2*H′1*QB′r′,
                                        ν^2*H1*Qr′]))
        else
            Qr = Qpref*H1
            QBr = QBpref*H1 
            Q′r = Qpref*H′1 
            QB′r = QBpref*H′1 
            #println("Q-prefactor = $Qpref")
            #println("Q′-prefactor = $QBpref")
            #println("Qr = $Qr, QBr = $QBr, Q′r = $Q′r, QB′r = $QB′r, H1r = $H1, H′1 = $H′1")
            #println("H1 = $H1")
            χ = SMatrix{3, 3}(im/(32*pi^2)*diagm([k^2*Q′r*H′1′ + ω^2*m^2/(r*r′*ν^2)*QBr*H1′,
                                        k^2*m^2/(ν^2*r*r′)*Qr*H1′ + ω^2*QB′r*H′1′,
                                        ν^2*Qr*H1′]))
        end
    end
    # Set NaN values to zero. So far these have only been seen for large im
    # ω such that response limits to zero.
    #χ[isnan.(χ)] .= 0
    # Set real response to real as soon as possible
    if real(ω) == 0
        return real.(χ)
    else 
        return χ
    end # if
end # function dynamicresponse
# Returns low-ν asymptote for the diff. response
function asymptdiffresponse(r, r′, ω, m, k)
    ν = getν(ω, k)
    if m == 0 #Not fixed yet
        return (32*pi^3)^-1*[k^2*r/r′; 0; 4im*k/r′;;
                             0; ω^2*r/r′; 0;;
                             0; 0; 0]
    else
        return (16*pi^3)^-1*(rw^2/(r*r′))^m *
        [m/(r*r′); -im*m/(r*r′); -im*k/r′;;
         im*m/(r*r′); m/(r*r′); k/r′;;
         im*k/r; k/r; ν^2/m]
    end
end

# Returns low-ν asymptote for the dynamic response
function asymptdynamicresponse(r, r′, ω, m, k)
    ν = getν(ω, k)
    if m == 0
        return (32*pi^3)^-1*[k^2*r/r′; 0; 4im*k/r′;;
                             0; ω^2*r/r′; 0;;
                             0; 0; 0]
    else
        if r < r′
            return (16*pi^3)^-1*(r/r′)^m*[m/(r*r′); im*m/(r*r′); im*k/r′;;
                                           im*m/(r*r′); -m/(r*r′); -k/r′;;
                                           im*k/r; -k/r; 0]
        else
            return (16*pi^3)^-1*(r′/r)^m* [m/(r*r′); -im*m/(r*r′); -im*k/r′;;
                                           -im*m/(r*r′); -m/(r*r′); -k/r′;;
                                           -im*k/r; -k/r; 0]
        end
    end
end

# Plots the local renormalized response vs ω for a given r
function localresponseωplot(r; ωmin=-20, ωstep=0.1, ωmax=20)
    ωω = ωmin:ωstep:ωmax
    χrr, χθr, χzr, χrθ, χθθ, χzθ, χrz, χθz, χzz = [[] for _ in 1:9]
    for ω in ωω

        χ = localresponse(r, ω + iϵ)
        push!(χrr, χ[1, 1])
        push!(χθr, χ[2, 1])
        push!(χzr, χ[3, 1])
        push!(χrθ, χ[1, 2])
        push!(χθθ, χ[2, 2])
        push!(χzθ, χ[3, 2])
        push!(χrz, χ[1, 3])
        push!(χθz, χ[2, 3])
        push!(χzz, χ[3, 3])
        println("ω = $ω has been integrated")
    end # for
    # Mask low values, for now
    ymin = 1e-10
    χrr[abs.(χrr) .< ymin] .= ymin
    χrθ[abs.(χrθ) .< ymin] .= ymin
    χrz[abs.(χrz) .< ymin] .= ymin
    χθr[abs.(χθr) .< ymin] .= ymin
    χθθ[abs.(χθθ) .< ymin] .= ymin
    χθz[abs.(χθz) .< ymin] .= ymin
    χzr[abs.(χzr) .< ymin] .= ymin
    χzθ[abs.(χzθ) .< ymin] .= ymin
    χzz[abs.(χzz) .< ymin] .= ymin
    χrrplot = plot(ωω, abs.(χrr), line_z=imag(log.(χrr)), 
                  label="χrr")
    χθrplot = plot(ωω, abs.(χθr), line_z=imag(log.(χθr)), 
                  label="χθr")
    χzrplot = plot(ωω, abs.(χzr), line_z=imag(log.(χzr)), 
                  label="χzr")
    χrθplot = plot(ωω, abs.(χrθ), line_z=imag(log.(χrθ)), 
                  label="χrθ")
    χθθplot = plot(ωω, abs.(χθθ), line_z=imag(log.(χθθ)), 
                  label="χθθ")
    χzθplot = plot(ωω, abs.(χzθ), line_z=imag(log.(χzθ)), 
                  label="χzθ")
    χrzplot = plot(ωω, abs.(χrz), line_z=imag(log.(χrz)), 
                  label="χrz")
    χθzplot = plot(ωω, abs.(χθz), line_z=imag(log.(χθz)), 
                  label="χθz")
    χzzplot = plot(ωω, abs.(χzz), line_z=imag(log.(χzz)), 
                  label="χzz")
    plot(χrrplot, χrθplot, χrzplot, χθrplot, χθθplot, χθzplot, χzrplot, χzθplot, χzzplot, layout=(3,3), title="Local response", yscale=:log10, size=(800, 500), clims=(0, pi))
end # function localresponseωplot

# Plots the local renormalized response vs imaginary ω for a given r
# Also plots a Drude functional form 1/iω * 1/(iω-τ-1)
function localresponseiωplot(r; ωmin=0, ωstep=0.1, ωmax=20, τ=1)
    ωω = im*(ωmin:ωstep:ωmax)
    χrr, χθr, χzr, χrθ, χθθ, χzθ, χrz, χθz, χzz = [[] for _ in 1:9]
    for ω in ωω

        χ = localresponse(r, ω + iϵ)
        push!(χrr, χ[1, 1])
        push!(χθr, χ[2, 1])
        push!(χzr, χ[3, 1])
        push!(χrθ, χ[1, 2])
        push!(χθθ, χ[2, 2])
        push!(χzθ, χ[3, 2])
        push!(χrz, χ[1, 3])
        push!(χθz, χ[2, 3])
        push!(χzz, χ[3, 3])
        #println("ω = $ω has been integrated")
    end # for
    # Mask low values, for now
    ymin = 1e-10
    χrr[abs.(χrr) .< ymin] .= ymin
    χrθ[abs.(χrθ) .< ymin] .= ymin
    χrz[abs.(χrz) .< ymin] .= ymin
    χθr[abs.(χθr) .< ymin] .= ymin
    χθθ[abs.(χθθ) .< ymin] .= ymin
    χθz[abs.(χθz) .< ymin] .= ymin
    χzr[abs.(χzr) .< ymin] .= ymin
    χzθ[abs.(χzθ) .< ymin] .= ymin
    χzz[abs.(χzz) .< ymin] .= ymin
    ωω = abs.(ωω) # Fixes axis for plotting
    
    # Computes drude form for comparison
    ceiling = maximum(abs.([χrr; χθθ; χzz]))
    println(ceiling)
    drude = map(x -> 1/x*1/(x-τ^-1), ωω)
    drude[drude .> ceiling] .= ceiling
    drude[drude .< ymin] .= ymin

    χrrplot = plot(ωω, abs.(χrr), line_z=imag(log.(χrr)), 
                  label="χrr")
    #plot!(ωω, drude, label="Drude form")
    χθrplot = plot(ωω, abs.(χθr), line_z=imag(log.(χθr)), 
                  label="χθr")
    χzrplot = plot(ωω, abs.(χzr), line_z=imag(log.(χzr)), 
                  label="χzr")
    χrθplot = plot(ωω, abs.(χrθ), line_z=imag(log.(χrθ)), 
                  label="χrθ")
    χθθplot = plot(ωω, abs.(χθθ), line_z=imag(log.(χθθ)), 
                  label="χθθ")
    χzθplot = plot(ωω, abs.(χzθ), line_z=imag(log.(χzθ)), 
                  label="χzθ")
    χrzplot = plot(ωω, abs.(χrz), line_z=imag(log.(χrz)), 
                  label="χrz")
    χθzplot = plot(ωω, abs.(χθz), line_z=imag(log.(χθz)), 
                  label="χθz")
    χzzplot = plot(ωω, abs.(χzz), line_z=imag(log.(χzz)), 
                  label="χzz")
    plot(χrrplot, χrθplot, χrzplot, χθrplot, χθθplot, χθzplot, χzrplot, χzθplot, χzzplot, layout=(3,3), title="Local response", yscale=:log10, xlabel="Im(ω)", size=(800, 500), clims=(0,pi))
end # function localresponseωplot

# Plots the ω integrand for drude response vs im. ω for a given r, r0
function drudeint1plot(r, r0; ymin=0, ystep=0.1, ymax=20, τ=1)
    yy = (ymin:ystep:ymax)
    inta, intb = [[] for _ in 1:2]

    χ0 = localresponse(r, iϵ)
    χτ = localresponse(r, im*τ^-1)

    a = √(r^2-r0^2)
    χ̃0 = a/r*χ0[1,1] + r0/(a*r)*χ0[2,2] .- r/a*χ0[3,3]
    χ̃τ = a/r*χτ[1,1] + rτ/(a*r)*χτ[2,2] .- r/a*χτ[3,3]

    for y in yy
        χ = localresponse(r, im*y)
        χ̃ = a/r*χ[1,1] + r0/(a*r)*χ[2,2] - r/a*χ[3,3]
        push!(inta, (χ̃ - χ̃0)/(y*(y-τ^-1)))
        push!(intb, (χ̃ - χ̃τ)/(y*(y-τ^-1)))
    end # for
    ymin = 1e-10

    # Mask low values, for now
    χrr[abs.(χrr) .< ymin] .= ymin
    χθθ[abs.(χθθ) .< ymin] .= ymin
    χzz[abs.(χzz) .< ymin] .= ymin
    ωω = abs.(ωω) # Fixes axis for plotting

    χrrplot = plot(ωω, abs.(χrr), line_z=imag(log.(χrr)), 
                  label="χrr")
    #plot!(ωω, drude, label="Drude form")
    χθrplot = plot(ωω, abs.(χθr), line_z=imag(log.(χθr)), 
                  label="χθr")
    χzrplot = plot(ωω, abs.(χzr), line_z=imag(log.(χzr)), 
                  label="χzr")
    χrθplot = plot(ωω, abs.(χrθ), line_z=imag(log.(χrθ)), 
                  label="χrθ")
    χθθplot = plot(ωω, abs.(χθθ), line_z=imag(log.(χθθ)), 
                  label="χθθ")
    χzθplot = plot(ωω, abs.(χzθ), line_z=imag(log.(χzθ)), 
                  label="χzθ")
    χrzplot = plot(ωω, abs.(χrz), line_z=imag(log.(χrz)), 
                  label="χrz")
    χθzplot = plot(ωω, abs.(χθz), line_z=imag(log.(χθz)), 
                  label="χθz")
    χzzplot = plot(ωω, abs.(χzz), line_z=imag(log.(χzz)), 
                  label="χzz")
    plot(χrrplot, χrθplot, χrzplot, χθrplot, χθθplot, χθzplot, χzrplot, χzθplot, χzzplot, layout=(3,3), title="Local response", yscale=:log10, xlabel="Im(ω)", size=(800, 500), clims=(0,pi))
end # function localresponseωplot

# Plots the dynamic, static and LSE response vs k for a given r
function kcheck(r, ω, m; kmin=0.0001, kstep=0.001, kmax=1, rd=1)
    kk = kmin:kstep:kmax
    #χsrr, χsθr, χszr, χsrθ, χsθθ, χszθ, χsrz, χsθz, χszz = [[] for _ in 1:9]
    χxrr, χxθr, χxzr, χxrθ, χxθθ, χxzθ, χxrz, χxθz, χxzz = [[] for _ in 1:9] # Dynamic
    χdrr, χdθr, χdzr, χdrθ, χdθθ, χdzθ, χdrz, χdθz, χdzz = [[] for _ in 1:9] # Diff.
    #Err, Eθr, Ezr, Erθ, Eθθ, Ezθ, Erz, Eθz, Ezz = [[] for _ in 1:9]
    #Brr, Bθr, Bzr, Brθ, Bθθ, Bzθ, Brz, Bθz, Bzz = [[] for _ in 1:9]
    for k in kk
#        α = dlayerresponsecoef(ω, m, k, [1; 0; 0], rd=rd)
#        LSEresponse = coeftoresponse(r, rd, ω, m, k, α)
#        push!(Err, LSEresponse[1])
#        push!(Eθr, LSEresponse[2])
#        push!(Ezr, LSEresponse[3])
#        push!(Brr, LSEresponse[4])
#        push!(Bθr, LSEresponse[5])
#        push!(Bzr, LSEresponse[6])
#
#        α = dlayerresponsecoef(ω, m, k, [0; 1; 0], rd=rd)
#        LSEresponse = coeftoresponse(r, rd, ω, m, k, α)
#        push!(Erθ, LSEresponse[1])
#        push!(Eθθ, LSEresponse[2])
#        push!(Ezθ, LSEresponse[3])
#        push!(Brθ, LSEresponse[4])
#        push!(Bθθ, LSEresponse[5])
#        push!(Bzθ, LSEresponse[6])
#
#        α = dlayerresponsecoef(ω, m, k, [0; 0; 1], rd=rd)
#        LSEresponse = coeftoresponse(r, rd, ω, m, k, α)
#        push!(Erz, LSEresponse[1])
#        push!(Eθz, LSEresponse[2])
#        push!(Ezz, LSEresponse[3])
#        push!(Brz, LSEresponse[4])
#        push!(Bθz, LSEresponse[5])
#        push!(Bzz, LSEresponse[6])

#        χs = staticresponse(r, rd, m, k)
#        push!(χsrr, χs[1, 1])
#        push!(χsθr, χs[2, 1])
#        push!(χszr, χs[3, 1])
#        push!(χsrθ, χs[1, 2])
#        push!(χsθθ, χs[2, 2])
#        push!(χszθ, χs[3, 2])
#        push!(χsrz, χs[1, 3])
#        push!(χsθz, χs[2, 3])
#        push!(χszz, χs[3, 3])

#        χx = dynamicresponse(r, rd, ω + iϵ, m, k)
#        push!(χxrr, χx[1, 1])
#        push!(χxθr, χx[2, 1])
#        push!(χxzr, χx[3, 1])
#        push!(χxrθ, χx[1, 2])
#        push!(χxθθ, χx[2, 2])
#        push!(χxzθ, χx[3, 2])
#        push!(χxrz, χx[1, 3])
#        push!(χxθz, χx[2, 3])
#        push!(χxzz, χx[3, 3])

        χd = diffresponse(r, rd, ω + iϵ, m, k)
        push!(χdrr, χd[1, 1])
        push!(χdθr, χd[2, 1])
        push!(χdzr, χd[3, 1])
        push!(χdrθ, χd[1, 2])
        push!(χdθθ, χd[2, 2])
        push!(χdzθ, χd[3, 2])
        push!(χdrz, χd[1, 3])
        push!(χdθz, χd[2, 3])
        push!(χdzz, χd[3, 3])
    end # for
    # Mask low values, for now
    ymin = 1e-8
    χdrr[abs.(χdrr) .< ymin] .= ymin
    χdrθ[abs.(χdrθ) .< ymin] .= ymin
    χdrz[abs.(χdrz) .< ymin] .= ymin
    χdθr[abs.(χdθr) .< ymin] .= ymin
    χdθθ[abs.(χdθθ) .< ymin] .= ymin
    χdθz[abs.(χdθz) .< ymin] .= ymin
    χdzr[abs.(χdzr) .< ymin] .= ymin
    χdzθ[abs.(χdzθ) .< ymin] .= ymin
    χdzz[abs.(χdzz) .< ymin] .= ymin
    χxrr[abs.(χxrr) .< ymin] .= ymin
    χxrθ[abs.(χxrθ) .< ymin] .= ymin
    χxrz[abs.(χxrz) .< ymin] .= ymin
    χxθr[abs.(χxθr) .< ymin] .= ymin
    χxθθ[abs.(χxθθ) .< ymin] .= ymin
    χxθz[abs.(χxθz) .< ymin] .= ymin
    χxzr[abs.(χxzr) .< ymin] .= ymin
    χxzθ[abs.(χxzθ) .< ymin] .= ymin
    χxzz[abs.(χxzz) .< ymin] .= ymin
#    χsrrplot = plot(kk, abs.(χsrr), #line_z=imag(log.(χsrr)), 
#                  label="χsrr")
#    χsθrplot = plot(kk, abs.(χsθr), #line_z=imag(log.(χsθr)), 
#                  label="χsθr")
#    χszrplot = plot(kk, abs.(χszr), #line_z=imag(log.(χszr)), 
#                  label="χszr")
#    χsrθplot = plot(kk, abs.(χsrθ), #line_z=imag(log.(χsrθ)), 
#                  label="χsrθ")
#    χsθθplot = plot(kk, abs.(χsθθ), #line_z=imag(log.(χsθθ)), 
#                  label="χsθθ")
#    χszθplot = plot(kk, abs.(χszθ), #line_z=imag(log.(χszθ)), 
#                  label="χszθ")
#    χsrzplot = plot(kk, abs.(χsrz), #line_z=imag(log.(χsrz)), 
#                  label="χsrz")
#    χsθzplot = plot(kk, abs.(χsθz), #line_z=imag(log.(χsθz)), 
#                  label="χsθz")
#    χszzplot = plot(kk, abs.(χszz), #line_z=imag(log.(χszz)), 
#                  label="χszz")
#    plot(χsrrplot, χsrθplot, χsrzplot, χsθrplot, χsθθplot, χsθzplot, χszrplot, χszθplot, χszzplot, layout=(3,3), title="Static", yscale=:log10)
#    χxrrplot = plot(kk, abs.(χxrr), line_z=imag(log.(χxrr)), 
#                  label="χxrr")
#    χxθrplot = plot(kk, abs.(χxθr), line_z=imag(log.(χxθr)), 
#                  label="χxθr")
#    χxzrplot = plot(kk, abs.(χxzr), line_z=imag(log.(χxzr)), 
#                  label="χxzr")
#    χxrθplot = plot(kk, abs.(χxrθ), line_z=imag(log.(χxrθ)), 
#                  label="χxrθ")
#    χxθθplot = plot(kk, abs.(χxθθ), line_z=imag(log.(χxθθ)), 
#                  label="χxθθ")
#    χxzθplot = plot(kk, abs.(χxzθ), line_z=imag(log.(χxzθ)), 
#                  label="χxzθ")
#    χxrzplot = plot(kk, abs.(χxrz), line_z=imag(log.(χxrz)), 
#                  label="χxrz")
#    χxθzplot = plot(kk, abs.(χxθz), line_z=imag(log.(χxθz)), 
#                  label="χxθz")
#    χxzzplot = plot(kk, abs.(χxzz), line_z=imag(log.(χxzz)), 
#                  label="χxzz")
#    plot(χxrrplot, χxrθplot, χxrzplot, χxθrplot, χxθθplot, χxθzplot, χxzrplot, χxzθplot, χxzzplot, layout=(3,3), title="Dynamic", yscale=:log10)
    χdrrplot = plot(kk, abs.(χdrr), line_z=imag(log.(χdrr)), 
                  label="χdrr")
    χdθrplot = plot(kk, abs.(χdθr), line_z=imag(log.(χdθr)), 
                  label="χdθr")
    χdzrplot = plot(kk, abs.(χdzr), line_z=imag(log.(χdzr)), 
                  label="χdzr")
    χdrθplot = plot(kk, abs.(χdrθ), line_z=imag(log.(χdrθ)), 
                  label="χdrθ")
    χdθθplot = plot(kk, abs.(χdθθ), line_z=imag(log.(χdθθ)), 
                  label="χdθθ")
    χdzθplot = plot(kk, abs.(χdzθ), line_z=imag(log.(χdzθ)), 
                  label="χdzθ")
    χdrzplot = plot(kk, abs.(χdrz), line_z=imag(log.(χdrz)), 
                  label="χdrz")
    χdθzplot = plot(kk, abs.(χdθz), line_z=imag(log.(χdθz)), 
                  label="χdθz")
    χdzzplot = plot(kk, abs.(χdzz), line_z=imag(log.(χdzz)), 
                  label="χdzz")
    plot(χdrrplot, χdrθplot, χdrzplot, χdθrplot, χdθθplot, χdθzplot, χdzrplot, χdzθplot, χdzzplot, layout=(3,3), title="Diff.", yscale=:log10)
#    Errplot = plot(kk, abs.(Err), #line_z=imag(log.(Err)), 
#                  label="Err")
#    Eθrplot = plot(kk, abs.(Eθr), #line_z=imag(log.(Eθr)), 
#                  label="Eθr")
#    Ezrplot = plot(kk, abs.(Ezr), #line_z=imag(log.(Ezr)), 
#                  label="Ezr")
#    Erθplot = plot(kk, abs.(Erθ), #line_z=imag(log.(Erθ)), 
#                  label="Erθ")
#    Eθθplot = plot(kk, abs.(Eθθ), #line_z=imag(log.(Eθθ)), 
#                  label="Eθθ")
#    Ezθplot = plot(kk, abs.(Ezθ), #line_z=imag(log.(Ezθ)), 
#                  label="Ezθ")
#    Erzplot = plot(kk, abs.(Erz), #line_z=imag(log.(Erz)), 
#                  label="Erz")
#    Eθzplot = plot(kk, abs.(Eθz), #line_z=imag(log.(Eθz)), 
#                  label="Eθz")
#    Ezzplot = plot(kk, abs.(Ezz), #line_z=imag(log.(Ezz)), 
#                  label="Ezz")
#    plot(Errplot, Erθplot, Erzplot, Eθrplot, Eθθplot, Eθzplot, Ezrplot, Ezθplot, Ezzplot, layout=(3,3), title="LSE E")
#    Brrplot = plot(kk, abs.(Brr), #line_z=imag(log.(Brr)), 
#                  label="Brr")
#    Bθrplot = plot(kk, abs.(Bθr), #line_z=imag(log.(Bθr)), 
#                  label="Bθr")
#    Bzrplot = plot(kk, abs.(Bzr), #line_z=imag(log.(Bzr)), 
#                  label="Bzr")
#    Brθplot = plot(kk, abs.(Brθ), #line_z=imag(log.(Brθ)), 
#                  label="Brθ")
#    Bθθplot = plot(kk, abs.(Bθθ), #line_z=imag(log.(Bθθ)), 
#                  label="Bθθ")
#    Bzθplot = plot(kk, abs.(Bzθ), #line_z=imag(log.(Bzθ)), 
#                  label="Bzθ")
#    Brzplot = plot(kk, abs.(Brz), #line_z=imag(log.(Brz)), 
#                  label="Brz")
#    Bθzplot = plot(kk, abs.(Bθz), #line_z=imag(log.(Bθz)), 
#                  label="Bθz")
#    Bzzplot = plot(kk, abs.(Bzz), #line_z=imag(log.(Bzz)), 
#                  label="Bzz")
    #plot(Brrplot, Brθplot, Brzplot, Bθrplot, Bθθplot, Bθzplot, Bzrplot, Bzθplot, Bzzplot, layout=(3,3), title="LSE B")
end # function kcheck

# Prints the discontinuities over rd and similar checks for the static and LSE
# solution
function disccheck(ω, m, k; rd=1, ϵ=0.00001)
    # LSE dz
    α = dlayerresponsecoef(ω, m, k, [0;0;1], rd=rd)
    LSEplus = coeftoresponse(rd + ϵ, rd, ω, m, k, α)
    LSEminus = coeftoresponse(rd - ϵ, rd, ω, m, k, α)
    println("The discontinuity over rd from the LSE (dz) is $(round.(LSEplus.-LSEminus, sigdigits=3))")
    # Static dz
    staticdiff = (staticresponse(rd+ϵ, rd, m, k)[:,3] .- staticresponse(rd-ϵ, rd, m, k)[:,3]) * (2*pi)^3
    println("The discontinuity over rd from the static solution is (dz) $(round.(staticdiff, sigdigits=3))")
    # Dynamic dz
    dyndiff = (dynamicresponse(rd+ϵ, rd, ω + iϵ, m, k)[:,3] .- dynamicresponse(rd-ϵ, rd, ω + iϵ, m, k)[:,3]) * (2*pi)^3
    println("The discontinuity over rd from the dynamic solution is (dz) $(round.(dyndiff, sigdigits=3))")
    # Real dz
    real_discz = 4pi*[-im*k/rd, 0, 0]
    println("This number should be $(round.(real_discz, sigdigits=3))")
    # dz at wire
    LSEw = coeftoresponse(rw, rd, ω, m, k, α)
    staticw = staticresponse(rw, rd, m, k)[:,3]
    println("The LSE (dz) gives at rw $LSEw")
    println("The static solution (dz) gives at rw $staticw")

    # LSE dθ
    α = dlayerresponsecoef(ω, m, k, [0;1;0], rd=rd)
    LSEplus = coeftoresponse(rd + ϵ, rd, ω, m, k, α)
    LSEminus = coeftoresponse(rd - ϵ, rd, ω, m, k, α)
    println("The discontinuity over rd from the LSE (dθ) is $(round.(LSEplus.-LSEminus, sigdigits=2))")
    # Static dθ
    staticdiff = (staticresponse(rd+ϵ, rd, m, k)[:,2] .- staticresponse(rd-ϵ, rd, m, k)[:,2])*(2*pi)^3
    println("The discontinuity over rd from the static solution (dθ) is $(round.(staticdiff, sigdigits=3))")
    # Dynamic dθ
    dyndiff = (dynamicresponse(rd+ϵ, rd, ω + iϵ, m, k)[:,2] .- dynamicresponse(rd-ϵ, rd, ω + iϵ, m, k)[:,2]) * (2*pi)^3
    println("The discontinuity over rd from the dynamic solution (dθ) is $(round.(dyndiff, sigdigits=3))")
    # Real dθ
    real_discθ = 4pi*[-im*m/rd^2, 0, 0]
    println("This number should be $(round.(real_discθ, sigdigits=3))")
    # dθ at wire
    LSEw = coeftoresponse(rw, rd, ω, m, k, α)
    staticw = staticresponse(rw, rd, m, k)[:,2]
    println("The LSE (dθ) gives at rw $LSEw")
    println("The static solution (dθ) gives at rw $staticw")

end # function disccheck

# Returns local response diff at position r
function localresponse(r, ω)
    maxk = 3*abs(ω) # This could maybe be refined for ω imaginary
    if real(ω) == 0
        maxk = 100*r^-1
    end # if
    #println(maxk)
    mcutq::Float64 = 1e-6
    mcuta::Float64 = 1e-50 # Absolute value m-cutoff to prevent runaway loop
    χ = @SMatrix zeros(Complex, 3, 3)
    χm = @SMatrix zeros(Complex, 3, 3)
    m::Int64 = 0
    while true
        #println("m = $m")
        χm = (quadgk((k -> diffresponse(r, r, ω + iϵ, m, k, diagonal=true)), -maxk, maxk)::Tuple{SMatrix{3, 3}, Float64})[1]
        # The above quadgk is type unstable, limiting performance

        χ::SMatrix{3, 3, Complex, 9} = χ + 2 *χm
        if m == 0
            χ = χ ./ 2
        end # if
        #println(χm)
        χmnorm::Float64 = norm(χm)
        if χmnorm/norm(χ) < mcutq || χmnorm < mcuta
            #println("Broke sum over m at m = $m")
            break
        end # if
        m += 1
    end # for
    return χ
end # function localresponse

# Integrates the A coefficient for ΔE due to the anis. part of an 
# anisotropic Drude ES response with relaxation time τ and conductivity 
# anis. δσ.
# The ES plane is taken to have smallest distance r0 to the wire centre
function ΔEdrudecoef(τ, δσ, r0)
    # Integration over x direction
    maxrq = 1e7
    rint = quadgk(r -> ΔEdrudeint1(r, τ, r0), r0, maxrq*r0)
    return -δσ/(32*pi*τ)*rint
end # function ΔEdrude

function ΔEdrudeint2(r, τ, r0)
    maxrq = 1e7
    intfull = quadgk(r -> ΔEdrudeint1(r, τ, r0), r0, maxrq*r0)[1]
    return intfull
end # function ΔEdrudeint1

function ΔEdrudeint1old(r, r0, τ)
    maxωq = 1e7 # Confirmed to be enough for (10, 1, 1)
    χ̃0 = ΔEdrudeχ̃(0, r, r0)
    χ̃τ = ΔEdrudeχ̃(im*τ^-1, r, r0)
    inta = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - χ̃0)/(y*(y-τ^-1)), -maxωq*τ^-1, 0, τ^-1/2)[1]
    println("First integral done")
    intb = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - χ̃τ)/(y*(y-τ^-1)), τ^-1/2, τ^-1,maxωq*τ^-1)[1]
    println("Second integral done")
    return inta + intb
end # function ΔEdrudeint1

function ΔEdrudeint1(r, r0, τ)
    α = 1000 # Lower cutoff for pole handling
    β = 10000 # Upper cutoff for pole handling
    maxωq = 1e6 # Confirmed to be enough for (10, 1, 1)


    χ̃0 = ΔEdrudeχ̃(0, r, r0)
    χ̃τ = ΔEdrudeχ̃(im*τ^-1, r, r0)
    χ̃0f(y) = y > -α*τ^-1 ? χ̃0 : 0
    χ̃τf(y) = y < β*τ^-1 ? χ̃τ : 0

    intb = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - χ̃τf(y))/(y*(y-τ^-1)), τ^-1/2, τ^-1,maxωq*τ^-1)[1]
    println("Second integral done")
    inta = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - χ̃0f(y))/(y*(y-τ^-1)), -maxωq*τ^-1, 0, τ^-1/2)[1]
    println("First integral done")
    vanish = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - χ̃τf(y))/(y*(y-τ^-1)), -1e8, -α*τ^-1)[1]
    println("Correction: $(χ̃τ*log(1 - β^-1) - χ̃0*log(1 + α^-1))")
    println("This should be vanishing: $vanish")
    return inta + intb + χ̃τ*log(1 - β^-1) - χ̃0*log(1 + α^-1)
end # function ΔEdrudeint1

function ΔEdrudeint(r, r0, τ)
    maxωq = 1e7
    maxrq = 1e7
    inta = hcubature(y, r -> (ΔEdrudeχ̃(im*y, r, r0) - ΔEdrudeχ̃(0, r, r0))/(y*(y-τ^-1)),
                     (-maxωq*τ^-1, r0), (0, τ^-1/2)[1])
    println("First integral done")
    intb = quadgk(y -> (ΔEdrudeχ̃(im*y, r, r0) - ΔEdrudeχ̃(im*τ^-1, r, r0))/(y*(y-τ^-1)), τ^-1/2, τ^-1,maxωq*τ^-1)[1]
    println("Second integral done")
    return inta + intb
end # function ΔEdrudeint1

# Integrands for ΔEdrudecoef, factor δσ/2 off from paper version
function ΔEdrudeχ̃(ω, r, r0)
    a = √(r^2-r0^2)
    χ = localresponse(r, ω)
    χ̃ = a/r*χ[1, 1] + r0/(a*r)*χ[2, 2] - r/a*χ[3, 3]
    return χ̃
end # function ΔEdrudeχ̃

# Plots integrands for ΔEdrudeint1 as a function of y
function ΔEdrudeχ̃plot(r, r0, τ; ymin=-5, ystep=0.01, ymax=5)
    yy = ymin:ystep:ymax
    χ̃0 = ΔEdrudeχ̃(0, r, r0)
    χ̃τ = ΔEdrudeχ̃(im*τ^-1, r, r0)
    integ1 = []
    integ2 = []

    println("χ̃0 = $χ̃0")
    println("χ̃τ = $χ̃τ")

    for y in yy
        χ̃ = ΔEdrudeχ̃(im*y, r, r0)

        push!(integ1, (χ̃-χ̃0)/(y*(y-τ^-1)))
        push!(integ2, (χ̃-χ̃τ)/(y*(y-τ^-1)))
        println("y = $y done!")
    end # for
    # Mask low values, for now
#    ymin = 1e-10
#    χzθ[abs.(χzθ) .< ymin] .= ymin
#    χzz[abs.(χzz) .< ymin] .= ymin

    integralplot = plot(yy, abs.(integ1), 
                        label="Integrand 1", size=(800, 500), title="Integrands")
    plot!(yy, abs.(integ2), 
                  label="Integrand 2", yscale=:log10)
    #plot!([τ^-1], seriestype=:vline, label="", color=:red, annotations=(τ^-1*1.1, minimum(abs.(integ1)), ("τ^-1'", 6, :bottom, :left, :red)))
end # function ΔEdrudeχ̃plot

# Integrates the energy shift due to an isotropic Drude ES response with relaxation time
# τ and base conductivity σ0. 
# The ES plane is taken to have smallest distance r0 to the wire centre
# This is according to the ´´old'' formula
function ΔEdrudeisoold(τ, σ0, r0)
    # Integrand for real space integration
    ωint = r -> (im*σ0*tr(localresponse(r, iϵ) - localresponse(r, im*τ^-1)))
    # Integration over non-z direction
    ΔE = quadgk(r -> r/√(r^2-r0^2)*ωint(r), r0, Inf)[1]

    return ΔE
end # function ΔEdrudeiso

# Integrates the A coefficient for ΔE due to the anis. part of an 
# anisotropic Drude ES response with relaxation time τ and conductivity 
# anis. δσ.
# The ES plane is taken to have smallest distance r0 to the wire centre
# This is according to the ´´old'' formula
function ΔEdrudecoefold(τ, δσ, r0)
    # Integration over x direction
    maxrq = 1e7
    A = imag(quadgk(r -> -im*δσ/8*ΔEdrudeintold(r, τ, r0), r0, r0*maxrq)[1]) # Currently factor i off
    return A
end # function ΔEdrude

# Integrands for ΔEdrudecoef, factor δσ/2 off from paper version
# This is according to the ´´old'' formula
function ΔEdrudeintold(r, τ, r0)
    a = √(r^2-r0^2)
    χ̃ = localresponse(r, iϵ) - localresponse(r, im*τ^-1)
    A = a/r*χ̃[1, 1] + r0^2/(r*a)*χ̃[2, 2] - r/a*χ̃[3, 3]
    return A
end # function ΔEdrudeint

# Plots (maximum) energy shift dependance due to the anis. part of Drude ES
# response as a function of relaxation time τ
# This is according to the ´´old'' formula
function drudeτplotold(δσ, r0; τmin=0.1, τstep=0.1, τmax=5)
    ττ = τmin:τstep:τmax
    EE = []
    for τ in ττ
        E = ΔEdrudecoefold(τ, δσ, r0)
        push!(EE, E)
        println("τ = $τ has been integrated!")
    end # for
    plot(ττ, EE, xlabel="τ", ylabel="Maximal energy shift due to wire direction", title="Energy shift due to local Drude χES, r0=$r0")
end # function drudeτplot

# Plots (maximum) energy shift dependance due to the anis. part of Drude ES
# response as a function of wire distance r0
# This is according to the ´´old'' formula
function druder0plotold(δσ, τ; r0min=0.1, r0step=0.1, r0max=5)
    rr = r0min:r0step:r0max
    EE = []
    for r0 in rr
        E = ΔEdrudecoefold(τ, δσ, r0)
        push!(EE, E)
        println("r0 = $r0 has been integrated!")
    end # for
    plot(rr, EE, xlabel="r0", ylabel="Maximal energy shift due to wire direction", title="Energy shift due to local Drude χES, τ=$τ")
end # function druder0plot

# Plots Hankel functions for clarity
function hankelplot(m; step=0.01, xmin=0, xmax=2)
    range = xmin:step:xmax
    H1 = hankelh1.(m, range)
    H2 = hankelh2.(m, range)
    H′1 = 1/2*(hankelh1.(m-1, range) - hankelh1.(m+1, range))
    H′2 = 1/2*(hankelh2.(m-1, range) - hankelh2.(m+1, range))

    plot(norm.(H1))
    plot!(norm.(H2))
    plot!(norm.(H′1))
    plot!(norm.(H′2))
end # function hankelplot

# Plots Q, QB, Qdiff
function qcheck(ω, m, k; minr=0.1, stepr=0.001, maxr=3)
    ν = getν(ω, k)
    rr = minr:stepr:maxr
    QQ, QB, QQd, QBd, JJ, QQdA, HH′1w, HH′2w, HH′q = [[] for _ in 1:9]
    for r in rr
        #Retrieve Hankel fcn values
        H1w = hankelh1(m, ν*rw)
        H2w = hankelh2(m, ν*rw)
        H′1w = 1/2*(hankelh1(m-1, ν*rw) - hankelh1(m+1, ν*rw))
        H′2w = 1/2*(hankelh2(m-1, ν*rw) - hankelh2(m+1, ν*rw))

        Jw = besselj(m, ν*rw)
        Yw = bessely(m, ν*rw)
        J′w = 1/2*(besselj(m-1, ν*rw) - besselj(m+1, ν*rw))
        Y′w = 1/2*(bessely(m-1, ν*rw) - bessely(m+1, ν*rw))

        H1 = hankelh1(m, ν*r)
        H2 = hankelh2(m, ν*r)
        Qr = H2 - H2w/H1w*H1
        QBr = H2 - H′2w/H′1w*H1
        Jx2 = 2*besselj(m, ν*r)
        QdAs = -2*besselj(m, ν*rw)/bessely(m, ν*rw)*bessely(m, ν*r)

        Qdiff = -2*(Jw/(Jw + im*Yw))*H1
        QBdiff = -2*(J′w/(J′w + im*Y′w))*H1

        push!(QQ, Qr)
        push!(QB, QBr)
        push!(JJ, Jx2)
        push!(QQd, Qdiff)
        push!(QBd, QBdiff)
        push!(QQdA, QdAs)c

        push!(HH′1w, H′1w)
        push!(HH′2w, H′2w)
        push!(HH′q, H′2w/H′1w)
    end # for

    plot(rr, [
              #abs.(QQ), 
              #abs.(QB), 
              #abs.(JJ), 
              #abs.(QQd), 
              abs.(QBd),
              abs.(QQdA),
              #abs.(HH′1w),
              #abs.(HH′1w),
              #abs.(HH′q)
              ], labels=[
                        #"QQ";; 
                        #"QB";;
                        #"JJ";;
                        #"QQd";;
                        "QBd";;
                        "QQdA";;
                        #"H′1w";;
                        #"H′2w";;
                        #"HH′q"
                        ],
              yscale=:log10
             )
end # function qcheck

# Plots asymptotic forms vs small ν
function asymptcheck(r, r′, m; minν=0.00000001, stepν=0.001, maxν=1)
    νν = minν:stepν:maxν
    νν = νν .+ iϵ
    HQ = []
    H = [0.0im]
    h = []
    Q = [0.0im]
    q = []
    H′ = []
    h′ = []
    H′num = []
    Q′ = []
    q′ = []
    Q′num = []
    hq = []
    H′Q = []
    h′q = []
    HQ′ = []
    hq′ = []
    H′Q′ = []
    h′q′ = []
    for ν in νν
        #Retrieve Hankel fcn values
        H1w = hankelh1(m, ν*rw)
        H2w = hankelh2(m, ν*rw)
        H1 = hankelh1(m, ν*r)
        H′1 = 1/2*(hankelh1(m-1, ν*r) - hankelh1(m+1, ν*r))
        H′1′ = 1/2*(hankelh1(m-1, ν*r′) - hankelh1(m+1, ν*r′))
        H′2′ = 1/2*(hankelh2(m-1, ν*r′) - hankelh2(m+1, ν*r′))
        H1′ = hankelh1(m, ν*r′)
        H2′ = hankelh2(m, ν*r′)
        Qr′ = H2′ - H2w/H1w*H1′
        Q′r′ = H′2′ - H2w/H1w*H′1′

        # Do numerical derivatives for debugging
        H′1num = (H1-last(H))/(stepν*r)
        Q′r′num = (Qr′-last(Q))/(stepν*r′)
        push!(H′num, H′1num)
        push!(Q′num, Q′r′num)

        push!(H, H1)
        push!(h, -im/pi*factorial(m-1)*(ν*r/2)^-m)
        push!(Q, Qr′)
        push!(q, 2/factorial(m)*(ν*r′/2)^m)
        push!(H′, H′1)
        push!(h′, im/(2*pi)*factorial(m)*(ν*r/2)^(-m-1))
        push!(Q′, Q′r′)
        push!(q′, 1/factorial(m-1)*(ν*r′/2)^(m-1))
        push!(HQ, H1*Qr′)
        push!(hq, -2im/(pi*m)*(r′/r)^m)
        push!(H′Q, H′1*Qr′)
        push!(h′q, 2*im/(ν*r*pi)*(r′/r)^m)
        push!(HQ′, H1*Q′r′)
        push!(hq′, -2*im/(ν*r′*pi)*(r′/r)^m)
        push!(H′Q′, H′1*Q′r′)
        push!(h′q′, 2*m*im/(pi*ν^2*r*r′)*(r′/r)^m)
    end

    # Fix for bootstrapping num. der.
    popfirst!(H)
    popfirst!(Q)
    H′num[1] = NaN
    Q′num[1] = NaN

    # Make ν:s real for plotting
    νν = real.(νν)

    Hplot = plot(νν, abs.(H), label="H", color=:red)
    plot!(νν, abs.(h), label="Asympt. h", color=:green, linestyle=:dash)

    Qplot = plot(νν, abs.(Q), label="Q", color=:red)
    plot!(νν, abs.(q), label="Asympt. q", color=:green, linestyle=:dash)
  
    phaseplot = plot(νν, imag.(log.(q./Q)), label="Rel. phase Q")
    H′plot = plot(νν, abs.(H′), label="H′", color=:red)
    plot!(νν, abs.(h′), label="Asympt. h′", color=:green, linestyle=:dash)
    #plot!(νν, abs.(H′num), label="num. h′", color=:violet, linestyle=:dot)

    Q′plot = plot(νν, abs.(Q′), label="Q′", color=:red)
    plot!(νν, abs.(q′), label="Asympt. q′", color=:green, linestyle=:dash)
    #plot!(νν, abs.(Q′num), label="num. q′", color=:violet, linestyle=:dot)

    HQplot = plot(νν, abs.(HQ), label="HQ", color=:red)
    plot!(νν, abs.(hq), label="Asympt. hq", color=:green, linestyle=:dash)

    H′Qplot = plot(νν, abs.(H′Q), label="H′Q", color=:red)
    plot!(νν, abs.(h′q), label="Asympt. h′q", color=:green, linestyle=:dash)

    HQ′plot = plot(νν, abs.(HQ′), label="HQ′", color=:red)
    plot!(νν, abs.(hq′), label="Asympt. hq′", color=:green, linestyle=:dash)

    H′Q′plot = plot(νν, abs.(H′Q′), label="H′Q′", color=:red)
    plot!(νν, abs.(h′q′), label="Asympt. h′q′", color=:green, linestyle=:dash)

    plot(Hplot, phaseplot, Qplot, H′plot, Q′plot, HQplot, H′Qplot, HQ′plot, H′Q′plot, layout=[5,4])
end


function getν(ω, k)
    k = abs(k) # Removes analyt. n k over k=0, but retains b.c.
    return im*√Complex(k-ω)*√Complex(ω+k)
end # function ν

# Wrapper for besselj returning zero for too large arguments
# Should only be paired with something going to 0 faster than besselj
function mybesselj(m, z)
    try 
        return besselj(m, z)
    catch e
        if e == SpecialFunctions.AmosException(2)
            return 0
        end
    end
end # function mybesselj

# Wrapper for bessely returning zero for too large arguments
# Should only be paired with something going to 0 faster than bessely
function mybessely(m, z)
    try 
        return bessely(m, z)
    catch e
        if e == SpecialFunctions.AmosException(2)
            return 0
        end
    end
end # function mybesselj

end # module MaterialLambShift
