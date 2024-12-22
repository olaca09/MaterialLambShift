# This file contains code for the birefringent plates calculations
# The initial goal is to replicate the results of Munday, Barash and Capasso
# (2005). Notation will coincide when possible.

# Fix for Unitful units used in Symbolics calculations
#include("UnifulSymbolicsFix.jl")

export ϵ2osc
export paperϵplot
export twoplateΩ 
export twoplateM
export Dplot
export Mplot
export Ddebug
export ∂lnD∂θfn
export ∂lnD∂θ

# Two-oscillator model parameters used in Munday, Barash and Capasso (2005)
const quartzII = (1.920, 1.350, 2.093e14u"Hz2π", 2.040e16u"Hz2π")
const quartz⊥ = (1.960, 1.377, 2.093e14u"Hz2π", 2.024e16u"Hz2π")
const calciteII = (5.300, 1.683, 2.691e14u"Hz2π", 1.660e16u"Hz2π")
const calcite⊥ = (6.300, 1.182, 2.691e14u"Hz2π", 2.134e16u"Hz2π")
const BaTiO3II = (3595, 4.128, 0.850e14u"Hz2π", 0.841e16u"Hz2π")
const BaTiO3⊥ = (145, 4.064, 0.850e14u"Hz2π", 0.896e16u"Hz2π")

export quartzII
export quartz⊥
export calciteII
export calcite⊥
export BaTiO3II
export BaTiO3⊥

"Returns the dielectric function at imaginary frequency iξ for an undamped two-oscillator model
with parameters par=(cIR, cUV, ωIR, ωUV)"
#function ϵ2osc(ξ::Unitful.Frequency, par::Tuple{Real, Real, Unitful.Frequency, Unitful.Frequency})
function ϵ2osc(ξ, par)
    cIR, cUV, ωIR, ωUV = par
    return 1 + cIR/(1 + (ξ/ωIR)^2) + cUV/(1 + (ξ/ωUV)^2)
end # function ϵ2osc

# Standard dielectric functions
ϵquartzII(ξ) = ϵ2osc(ξ, quartzII)
ϵquartz⊥(ξ) = ϵ2osc(ξ, quartz⊥)
ϵcalciteII(ξ) = ϵ2osc(ξ, calciteII)
ϵcalcite⊥(ξ) = ϵ2osc(ξ, calcite⊥)
ϵBaTiO3II(ξ) = ϵ2osc(ξ, BaTiO3II)
ϵBaTiO3⊥(ξ) = ϵ2osc(ξ, BaTiO3⊥)
ϵ0(ξ) = 1

export ϵquartzII
export ϵquartz⊥
export ϵcalciteII
export ϵcalcite⊥
export ϵBaTiO3II
export ϵBaTiO3⊥
export ϵ0

"Plots the dielectric function for imaginary frequencies iξ in the log range 10^(ξmin, ξstep, ξmax), for
two-oscillator models with parameters as in the paper of Munday, Barash and Capasso (2005)"
function paperϵplot(ξmin, ξstep, ξmax)
    ξrange = 10 .^(ξmin:ξstep:ξmax)
    quartplot = plot(ξrange, ϵ2osc.(ξrange, (quartzII,)), label="Quartz∥", yscale=:identity)
    plot!(ξrange, ϵ2osc.(ξrange, (quartz⊥,)), label="Quartz⊥", yscale=:identity)
    calciteplot = plot(ξrange, ϵ2osc.(ξrange, (calciteII,)), label="Calcite∥", yscale=:identity)
    plot!(ξrange, ϵ2osc.(ξrange, (calcite⊥,)), label="Calcite⊥", yscale=:identity)
    BaTiO3plot = plot(ξrange, ϵ2osc.(ξrange, (BaTiO3II,)), label="BaTiO3∥", yscale=:log10)
    plot!(ξrange, ϵ2osc.(ξrange, (BaTiO3⊥,)), label="BaTiO3⊥", yscale=:log10)

    plot(quartplot, calciteplot, BaTiO3plot, layout=(3,1), size=(800, 600), xscale=:log10, xlabel="ξ", ylabel="ϵ(iξ)")
    savefig("extra_plots/paperϵplot.svg")
end # function papperϵplot

#"Returns the torque per unit area of the two-plate system with dielectric functions 
#ϵ1II(ξ), ϵ1⊥(ξ), ϵ2II(ξ) and ϵ2⊥(ξ), ϵ3(ξ), separation dval, principal axes angle θval and temperature T"
#function twoplateM(dval, θval, T; ξcutoff = 1e18u"Hz2π", rcutoff = 1e10u"m^-1")
#    # r is the transverse wavenumber
#
#    if !(typeof(T) <: Unitful.Temperature)
#        T = T*u"K"
#        println("Temperature was not given in temperature Unitful type, assuming Kelvin")
#    end # if

    @variables θ, r, ϕ, ξ, par1II[1:4], par1⊥[1:4], par2II[1:4], par2⊥[1:4], d, c
    export θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c


    # Definitions of integrand and parts thereof

    # The dielectric functions at ξ
    ϵ1IIξ = ϵ2osc(ξ, par1II)
    ϵ1⊥ξ = ϵ2osc(ξ, par1⊥)
    ϵ2IIξ = ϵ2osc(ξ, par2II)
    ϵ2⊥ξ = ϵ2osc(ξ, par2⊥)
    ϵ3ξ = ϵ0(ξ)

    # Store trig values for ϕ
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # Wave numbers in the z direction
    ρ1 = (r^2 + c^-2*ξ^2*ϵ1⊥ξ)^(1/2)
    ρ2 = (r^2 + c^-2*ξ^2*ϵ2⊥ξ)^(1/2)
    ρ3 = (r^2 + c^-2*ξ^2*ϵ3ξ)^(1/2)
    ρ̃1 = (r^2 + (ϵ1IIξ/ϵ1⊥ξ - 1)r^2 * cosϕ^2 + c^-2*ξ^2*ϵ1IIξ)^(1/2)
    ρ̃2 = (r^2 + (ϵ2IIξ/ϵ2⊥ξ - 1)*(r*cosϕ*cos(θ) - r*sinϕ*sin(θ))^2
          + c^-2*ξ^2*ϵ2IIξ)^(1/2)

    # Gruesome algebraic expressions
    γ = (
         (ρ1+ρ3) * (ρ2+ρ3) * (
                             (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) - ϵ1⊥ξ*(ρ̃1-ρ1)*(r^2*sinϕ^2-ρ1*ρ3)
                             /(ρ1^2-r^2*sinϕ^2)
                            )
        *(
          (ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - ϵ2⊥ξ*(ρ̃2-ρ2) * ((r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2-ρ2*ρ3)
          /(ρ2^2 - (r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2)
         )
       )

    A = (
         (
         (ρ1+ρ3)*(ρ2+ρ3) - (ρ1-ρ3)*(ρ2-ρ3)*exp(-2ρ3*d)
         )
        *(
          (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)*exp(-2ρ3*d)
        )
        - (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
        *(
          (r^2*sinϕ^2 - ρ1*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) * (ρ2+ρ3) * (ρ1+ρ3)
          + 2*(ϵ2⊥ξ-ϵ3ξ)*(
                          r^2*sinϕ^2 * (r^2*ρ1-ρ2*ρ3^2) 
                          + ρ1*ρ3^2*(r^2 - 2*r^2*sinϕ^2 + ρ1*ρ2)
                         )
          * exp(-2ρ3*d) + (r^2*sinϕ^2+ρ1*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)
          * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4ρ3*d)
        )
       )
    
    B = (
         (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*(ϵ1⊥ξ-ϵ3ξ) * (r^2*ρ2-ρ1*ρ3^2-2*ρ2*ρ3^2)
         * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
        )
        + (ρ̃1+ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
        *(
          -1*(r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3)*(ρ2+ρ3) + 2*(
                                                         r^2*sinϕ^2 *(ρ1*ρ2+ρ3^2)
                                                         - ρ1^2*ρ3^2 + ρ1*ρ2*ρ3^2
                                                        ) * exp(-2*ρ3*d)
          - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
         )

    C = (ρ2*ρ3*(
               -1*(ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*ρ3*(ϵ1⊥ξ-ϵ3ξ)
               * (r^2+ρ1*ρ2) * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3)
               * (ρ2-ρ3) * exp(-4*ρ3*d)
              )
    + (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)*ρ2*ρ3
    *(
      (r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3) * (ρ2+ρ3)
      + 2*ρ3*(
              ρ1^2*ρ2 + ρ1*ρ3^2 + r^2*sinϕ^2 * (ρ1-ρ2)
             ) * exp(-2*ρ3*d)
      - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
     )
   )

    E = 4*ρ1*ρ2*ρ3^2 * (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2) * exp(-2*ρ3*d)

    D = (γ^-1 * (A - (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
                        *(
                          B*r^2*sin(ϕ+θ)^2 - E*(
                                                2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                               )
                          + C
                         )
                       )
           )
    lnD = log(D)

    # Define derivative of lnD with respect to θ
    ∂θ = Differential(θ)
    ∂lnD∂θ = expand_derivatives(∂θ(lnD))

    # Build Julia function for ∂lnD∂θ
#    function ∂lnD∂θfn(θval, rval, ϕval, ξval, par1IIval, par1⊥val, par2IIval, par2⊥val, dval, cval)
#        return build_function(∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
#                              expression = Val{false})
#    end # function ∂lnD∂θfn
    ∂lnD∂θfn = eval(build_function(∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                                   expression = Val{true}))


    #println(expand_derivatives(∂lnD∂θ))

#    ## Try evaluating ∂lnD∂θ at some point
#    valuedict = Dict([θ => θval, r => 50119u"m^-1", ϕ => 0.58, ξ => 0u"Hz2π", par1II => calciteII, 
#                      par1⊥ => calcite⊥, par2II => BaTiO3II, par2⊥ => BaTiO3⊥, 
#                      d => dval, c => c_0])
#
#    println("∂lnD∂θ = $(substitute(∂lnD∂θ, valuedict))")
#
"Returns the torque per unit area of the two-plate system with dielectric functions 
ϵ1II(ξ), ϵ1⊥(ξ), ϵ2II(ξ) and ϵ2⊥(ξ), ϵ3(ξ), separation dval, principal axes angle θval and temperature T"
function twoplateM(dval, θval, T; ξcutoff = 1e18u"Hz2π", rcutoff = 1e10u"m^-1", rtol = 1e-2)
    # r is the transverse wavenumber

    if !(typeof(T) <: Unitful.Temperature)
        T = T*u"K"
        println("Temperature was not given in temperature Unitful type, assuming Kelvin")
    end # if

    # Matsubara sum #
    M = 0
    Merror = 0
    ξn = 0u"Hz2π"

    while ξn < ξcutoff
#        integ = rϕ -> rϕ[1] * substitute(∂lnD∂θ, Dict([θ => θval, r => rϕ[1]*u"m^-1", ϕ => rϕ[2], 
#                                                    ξ => ξn, par1II => calciteII, 
#                                                    par1⊥ => calcite⊥, par2II => BaTiO3II,
#                                                    par2⊥ => BaTiO3⊥, d => dval, c => c_0]))
    function integ(rϕ)
        #println("Calculating torque for r = $(rϕ[1]), ξ=$ξn, and ϕ = $(rϕ[2])")
        return rϕ[1] * ∂lnD∂θfn(θval, rϕ[1]*u"m^-1", rϕ[2], ξn, calciteII, calcite⊥,
                                        BaTiO3II, BaTiO3⊥, dval, c_0)
    end # function integ
        rcutoff = ustrip(u"m^-1", rcutoff) # Strip units for hcubature

        # Split integration into two domains in ϕ, to avoid divisions by zero
        # at sin(ϕ) = ± 1
        M, Merror += k_B*T/(4*π^2) * hcubature(integ, (0, -π/2), (rcutoff, π/2), rtol=rtol)*u"m^-2"
        M, Merror += k_B*T/(4*π^2) * hcubature(integ, (0, π/2), (rcutoff, 3π/2), rtol=rtol)*u"m^-2"
        # Unit of r dr is missing

        if ξn == 0 # Halve first term
            M = M/2
        end # if

        println("M contribution from ξn = $ξn calculated, current torque M = $M")
        ξn += uconvert(u"Hz2π", k_B*T/ħ)
    end # while

end # function twoplateΩ

"Returns the free energy per unit area of the two-plate system with dielectric functions 
ϵ1II(ξ), ϵ1⊥(ξ), ϵ2II(ξ) and ϵ2⊥(ξ), ϵ3(ξ), separation d, principal axes angle θ and temperature T"
function twoplateΩ(ϵ1II, ϵ1⊥, ϵ2II, ϵ2⊥, ϵ3, d, θ, T; ξcutoff = 1e18u"Hz2π", rcutoff = 1e10u"m^-1")
    # r is the transverse wavenumbero

    if !(typeof(T) <: Unitful.Temperature)
        T = T*u"K"
        println("Temperature was not given in temperature Unitful type, assuming Kelvin")
    end # if

    # Definitions of integrand and parts thereof
    function lnD(rϕ, ξ)

        r, ϕ = rϕ

        # Store the dielectric functions at ξ
        ϵ1IIξ = ϵ1II(ξ)
        ϵ1⊥ξ = ϵ1⊥(ξ)
        ϵ2IIξ = ϵ2II(ξ)
        ϵ2⊥ξ = ϵ2⊥(ξ)
        ϵ3ξ = ϵ3(ξ)

        # Store trig values for ϕ
        sinϕ = sin(ϕ)
        cosϕ = cos(ϕ)

        # Wave numbers in the z direction
        ρ1 = √(r^2 + c_0^-2*ξ^2*ϵ1⊥ξ)
        ρ2 = √(r^2 + c_0^-2*ξ^2*ϵ2⊥ξ)
        ρ3 = √(r^2 + c_0^-2*ξ^2*ϵ3ξ)
        ρ̃1 = √(r^2 + (ϵ1IIξ/ϵ1⊥ξ - 1)r^2 * cosϕ^2 + c_0^-2*ξ^2*ϵ1IIξ)
        ρ̃2 = √(r^2 + (ϵ2IIξ/ϵ2⊥ξ - 1)*(r*cosϕ*cos(θ) - r*sinϕ*sin(θ))^2
                                           + c_0^-2*ξ^2*ϵ2IIξ)

        # Gruesome algebraic expressions
        γ = (
             (ρ1+ρ3) * (ρ2+ρ3) * (
                                 (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) - ϵ1⊥ξ*(ρ̃1-ρ1)*(r^2*sinϕ^2-ρ1*ρ3)
                                 /(ρ1^2-r^2*sinϕ^2)
                                )
            *(
              (ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - ϵ2⊥ξ*(ρ̃2-ρ2) * ((r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2-ρ2*ρ3)
              /(ρ2^2 - (r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2)
             )
           )

        A = (
             (
             (ρ1+ρ3)*(ρ2+ρ3) - (ρ1-ρ3)*(ρ2-ρ3)*exp(-2ρ3*d)
             )
            *(
              (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)*exp(-2ρ3*d)
            )
            - (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
            *(
              (r^2*sinϕ^2 - ρ1*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) * (ρ2+ρ3) * (ρ1+ρ3)
              + 2*(ϵ2⊥ξ-ϵ3ξ)*(
                              r^2*sinϕ^2 * (r^2*ρ1-ρ2*ρ3^2) 
                              + ρ1*ρ3^2*(r^2 - 2*r^2*sinϕ^2 + ρ1*ρ2)
                             )
              * exp(-2ρ3*d) + (r^2*sinϕ^2+ρ1*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)
              * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4ρ3*d)
            )
           )
        
        B = (
             (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*(ϵ1⊥ξ-ϵ3ξ) * (r^2*ρ2-ρ1*ρ3^2-2*ρ2*ρ3^2)
             * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
            )
            + (ρ̃1+ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
            *(
              -1*(r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3)*(ρ2+ρ3) + 2*(
                                                             r^2*sinϕ^2 *(ρ1*ρ2+ρ3^2)
                                                             - ρ1^2*ρ3^2 + ρ1*ρ2*ρ3^2
                                                            ) * exp(-2*ρ3*d)
              - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
             )

        C = (ρ2*ρ3*(
                   -1*(ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*ρ3*(ϵ1⊥ξ-ϵ3ξ)
                   * (r^2+ρ1*ρ2) * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3)
                   * (ρ2-ρ3) * exp(-4*ρ3*d)
                  )
        + (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)*ρ2*ρ3
        *(
          (r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3) * (ρ2+ρ3)
          + 2*ρ3*(
                  ρ1^2*ρ2 + ρ1*ρ3^2 + r^2*sinϕ^2 * (ρ1-ρ2)
                 ) * exp(-2*ρ3*d)
          - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
         )
       )

        E = 4*ρ1*ρ2*ρ3^2 * (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2) * exp(-2*ρ3*d)
        logarg = ( γ^-1 * (A - (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
                            *(
                              B*r^2*sin(ϕ+θ)^2 - E*(
                                                    2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                                   )
                              + C
                             )
                           )
               )
        if logarg < 0
            println("D = $logarg")
            println("A = $A")
            println("B = $B")
            println("C = $C")
            println("E = $E")
            println("γ = $γ")
            println("ξ = $ξ")
            println("r = $r")
            println("ϕ = $ϕ")
            println("ρ1 = $ρ1")
            println("ρ2 = $ρ2")
            println("ρ3 = $ρ3")
            println("ρ̃1 = $ρ̃1")
            println("ρ̃2 = $ρ̃2")
        end # if


        return log( γ^-1 * (A - (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
                            *(
                              B*r^2*sin(ϕ+θ)^2 - E*(
                                                    2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                                   )
                              + C
                             )
                           )
                  )
    end # function lnD

    # Matsubara sum #
    Ω = 0
    Ωerror = 0
    ξn = 0u"Hz2π"

    while ξn < ξcutoff
        integ = rϕ -> rϕ[1] * lnD((rϕ[1]*u"m^-1", rϕ[2]), ξn)
        rcutoff = ustrip(u"m^-1", rcutoff) # Strip units for hcubature

        # Split integration into two domains in ϕ, to avoid divisions by zero
        # at sin(ϕ) = ± 1
        Ω, Ωerror += k_B*T/(4*π^2) * hcubature(integ, (0, -π/2), (rcutoff, π/2)) *u"m^-2" 
        Ω, Ωerror += k_B*T/(4*π^2) * hcubature(integ, (0, π/2), (rcutoff, 3π/2)) *u"m^-2" 
        # Unit of r dr is missing

        if ξn == 0 # Halve first term
            Ω = Ω/2
        end # if

        ξn += uconvert(u"Hz2π", k_B*T/ħ)
    end # while

end # function twoplateΩ

"Plots the torque contribution for twoplateM, over r and ϕ for a given ξ"
function Mplot(ξ, par1II, par1⊥, par2II, par2⊥, dval, θval; rmin = -4, rmax=15, ϕmin=0, ϕmax=2π, ϕstep=0.1, rstep=0.5)
    rrange = (10 .^ (rmin:rstep:rmax))u"m^-1"
    ϕrange = ϕmin:ϕstep:ϕmax
    MM = zeros(length(rrange), length(ϕrange))
    for i in 1:length(rrange)
        for j in 1:length(ϕrange)
            MM[i, j] = ∂lnD∂θfn(θval, rrange[i], ϕrange[j], ξ, par1II, par1⊥, par2II, par2⊥, dval, c_0)
        end # for
    end # for

    heatmap(ϕrange, rrange, MM, xlabel="ϕ", ylabel="r", title="∂lnD∂θ(r, ϕ) for ξ = $ξ",
           yscale=:log10)
end # function Mplot
"Plots the argument to the logarithm for twoplateΩ, D, over r and ϕ for a given ξ"
function Dplot(ξ, ϵ1II, ϵ1⊥, ϵ2II, ϵ2⊥, ϵ3, d, θ; rmin = -4, rmax=15, ϕmin=0, ϕmax=2π)
    rrange = (10 .^ (rmin:0.1:rmax))u"m^-1"
    ϕrange = ϕmin:0.01:ϕmax
    DD = zeros(length(rrange), length(ϕrange))

    # Definitions of integrand and its parts
    function D(r, ϕ, ξ)

        #println("Calculating D for r = $r and ϕ = $ϕ")

        # Store the dielectric functions at ξ
        ϵ1IIξ = ϵ1II(ξ)
        ϵ1⊥ξ = ϵ1⊥(ξ)
        ϵ2IIξ = ϵ2II(ξ)
        ϵ2⊥ξ = ϵ2⊥(ξ)
        ϵ3ξ = ϵ3(ξ)

        # Store trig values for ϕ
        sinϕ = sin(ϕ)
        cosϕ = cos(ϕ)

        # Wave numbers in the z direction
        ρ1 = √(r^2 + c_0^-2*ξ^2*ϵ1⊥ξ)
        ρ2 = √(r^2 + c_0^-2*ξ^2*ϵ2⊥ξ)
        ρ3 = √(r^2 + c_0^-2*ξ^2*ϵ3ξ)
        ρ̃1 = √(r^2 + (ϵ1IIξ/ϵ1⊥ξ - 1)r^2 * cosϕ^2 + c_0^-2*ξ^2*ϵ1IIξ)
        ρ̃2 = √(r^2 + (ϵ2IIξ/ϵ2⊥ξ - 1)*r^2*cos(ϕ+θ)^2
                                           + c_0^-2*ξ^2*ϵ2IIξ)

        # Gruesome algebraic expressions
        γ = (
             (ρ1+ρ3) * (ρ2+ρ3) * (
                                 (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) - ϵ1⊥ξ*(ρ̃1-ρ1)*(r^2*sinϕ^2-ρ1*ρ3)
                                 /(ρ1^2-r^2*sinϕ^2)
                                )
            *(
              (ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - ϵ2⊥ξ*(ρ̃2-ρ2) * ((r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2-ρ2*ρ3)
              /(ρ2^2 - (r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2)
             )
           )

        A = (
             (
             (ρ1+ρ3)*(ρ2+ρ3) - (ρ1-ρ3)*(ρ2-ρ3)*exp(-2ρ3*d)
             )
            *(
              (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)*exp(-2ρ3*d)
            )
            - (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
            *(
              (r^2*sinϕ^2 - ρ1*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) * (ρ2+ρ3) * (ρ1+ρ3)
              + 2*(ϵ2⊥ξ-ϵ3ξ)*(
                              r^2*sinϕ^2 * (r^2*ρ1-ρ2*ρ3^2) 
                              + ρ1*ρ3^2*(r^2 - 2*r^2*sinϕ^2 + ρ1*ρ2)
                             )
              * exp(-2ρ3*d) + (r^2*sinϕ^2+ρ1*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)
              * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4ρ3*d)
            )
           )
        
        B = (
             (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*(ϵ1⊥ξ-ϵ3ξ) * (r^2*ρ2-ρ1*ρ3^2-2*ρ2*ρ3^2)
             * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
            )
            + (ρ̃1+ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
            *(
              -1*(r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3)*(ρ2+ρ3) + 2*(
                                                             r^2*sinϕ^2 *(ρ1*ρ2+ρ3^2)
                                                             - ρ1^2*ρ3^2 + ρ1*ρ2*ρ3^2
                                                            ) * exp(-2*ρ3*d)
              - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
             )

        C = (ρ2*ρ3*(
                   -1*(ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*ρ3*(ϵ1⊥ξ-ϵ3ξ)
                   * (r^2+ρ1*ρ2) * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3)
                   * (ρ2-ρ3) * exp(-4*ρ3*d)
                  )
        + (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)*ρ2*ρ3
        *(
          (r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3) * (ρ2+ρ3)
          + 2*ρ3*(
                  ρ1^2*ρ2 + ρ1*ρ3^2 + r^2*sinϕ^2 * (ρ1-ρ2)
                 ) * exp(-2*ρ3*d)
          - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
         )
       )

        E = 4*ρ1*ρ2*ρ3^2 * (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2) * exp(-2*ρ3*d)
        logarg = ( γ^-1 * (A - (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
                            *(
                              B*r^2*sin(ϕ+θ)^2 - E*(
                                                    2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                                   )
                              + C
                             )
                           )
               )
        prefac_debug = (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
        curlybracket_debug = (B*r^2*sin(ϕ+θ)^2 - E*(
                                                    2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                                   )
                              + C)
        if logarg < 0
            println("")
            println("D = $logarg")
            println("A = $A")
            println("B = $B")
            println("C = $C")
            println("E = $E")
            println("γ = $γ")
            println("curly bracket = $curlybracket_debug")
            println("prefac = $prefac_debug")
            #println("ξ = $ξ")
            println("r = $r")
            println("ϕ = $ϕ")
            println("ρ1 = $ρ1")
            println("ρ2 = $ρ2")
            println("ρ3 = $ρ3")
            println("ρ̃1 = $ρ̃1")
            println("ρ̃2 = $ρ̃2")
        end # if

        return logarg
    end # function D
    for i in 1:length(rrange)
        for j in 1:length(ϕrange)
        #println("Called D for r = $(rrange[i]) and ϕ = $(ϕrange[j])")
            DD[i, j] = D(rrange[i], ϕrange[j], ξ)
        end # for
    end # for

    heatmap(ϕrange, rrange, DD, xlabel="ϕ", ylabel="r", title="D(r, ϕ) for ξ = $ξ",
           yscale=:log10)

#plot(ϕrange, DD[end, :], label="r = $(rrange[1])")
#plot(rrange, DD[:, end], label="ϕ = $(ϕrange[end])")
end # function Dplot

"Debug D"
function Ddebug(r, ϕ, ξ, ϵ1II, ϵ1⊥, ϵ2II, ϵ2⊥, ϵ3, d, θ)
    #println("Calculating D for r = $r and ϕ = $ϕ")

    # Store the dielectric functions at ξ
    ϵ1IIξ = ϵ1II(ξ)
    ϵ1⊥ξ = ϵ1⊥(ξ)
    ϵ2IIξ = ϵ2II(ξ)
    ϵ2⊥ξ = ϵ2⊥(ξ)
    ϵ3ξ = ϵ3(ξ)

    # Store trig values for ϕ
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # Wave numbers in the z direction
    ρ1 = √(r^2 + c_0^-2*ξ^2*ϵ1⊥ξ)
    ρ2 = √(r^2 + c_0^-2*ξ^2*ϵ2⊥ξ)
    ρ3 = √(r^2 + c_0^-2*ξ^2*ϵ3ξ)
    ρ̃1 = √(r^2 + (ϵ1IIξ/ϵ1⊥ξ - 1)r^2 * cosϕ^2 + c_0^-2*ξ^2*ϵ1IIξ)
    ρ̃2 = √(r^2 + (ϵ2IIξ/ϵ2⊥ξ - 1)*(r*cosϕ*cos(θ) - r*sinϕ*sin(θ))^2
                                       + c_0^-2*ξ^2*ϵ2IIξ)

    # Gruesome algebraic expressions
    γ = (
         (ρ1+ρ3) * (ρ2+ρ3) * (
                             (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) - ϵ1⊥ξ*(ρ̃1-ρ1)*(r^2*sinϕ^2-ρ1*ρ3)
                             /(ρ1^2-r^2*sinϕ^2)
                            )
        *(
          (ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - ϵ2⊥ξ*(ρ̃2-ρ2) * ((r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2-ρ2*ρ3)
          /(ρ2^2 - (r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2)
         )
       )

    A = (
         (
         (ρ1+ρ3)*(ρ2+ρ3) - (ρ1-ρ3)*(ρ2-ρ3)*exp(-2ρ3*d)
         )
        *(
          (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) - (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)*exp(-2ρ3*d)
        )
        - (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
        *(
          (r^2*sinϕ^2 - ρ1*ρ3)*(ϵ3ξ*ρ2+ϵ2⊥ξ*ρ3) * (ρ2+ρ3) * (ρ1+ρ3)
          + 2*(ϵ2⊥ξ-ϵ3ξ)*(
                          r^2*sinϕ^2 * (r^2*ρ1-ρ2*ρ3^2) 
                          + ρ1*ρ3^2*(r^2 - 2*r^2*sinϕ^2 + ρ1*ρ2)
                         )
          * exp(-2ρ3*d) + (r^2*sinϕ^2+ρ1*ρ3)*(ϵ3ξ*ρ2-ϵ2⊥ξ*ρ3)
          * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4ρ3*d)
        )
       )
    
    B = (
         (ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*(ϵ1⊥ξ-ϵ3ξ) * (r^2*ρ2-ρ1*ρ3^2-2*ρ2*ρ3^2)
         * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
        )
        + (ρ̃1+ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)
        *(
          -1*(r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3)*(ρ2+ρ3) + 2*(
                                                         r^2*sinϕ^2 *(ρ1*ρ2+ρ3^2)
                                                         - ρ1^2*ρ3^2 + ρ1*ρ2*ρ3^2
                                                        ) * exp(-2*ρ3*d)
          - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
         )

    C = (ρ2*ρ3*(
               -1*(ϵ3ξ*ρ1+ϵ1⊥ξ*ρ3) * (ρ1+ρ3) * (ρ2+ρ3) + 2*ρ3*(ϵ1⊥ξ-ϵ3ξ)
               * (r^2+ρ1*ρ2) * exp(-2*ρ3*d) + (ϵ3ξ*ρ1-ϵ1⊥ξ*ρ3) * (ρ1-ρ3)
               * (ρ2-ρ3) * exp(-4*ρ3*d)
              )
    + (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2)*ρ2*ρ3
    *(
      (r^2*sinϕ^2-ρ1*ρ3) * (ρ1+ρ3) * (ρ2+ρ3)
      + 2*ρ3*(
              ρ1^2*ρ2 + ρ1*ρ3^2 + r^2*sinϕ^2 * (ρ1-ρ2)
             ) * exp(-2*ρ3*d)
      - (r^2*sinϕ^2+ρ1*ρ3) * (ρ1-ρ3) * (ρ2-ρ3) * exp(-4*ρ3*d)
     )
   )

    E = 4*ρ1*ρ2*ρ3^2 * (ρ̃1-ρ1)*ϵ1⊥ξ/(ρ1^2-r^2*sinϕ^2) * exp(-2*ρ3*d)
    logarg = ( γ^-1 * (A - (ρ̃2-ρ2)*ϵ2⊥ξ/(ρ2^2-r^2*sin(ϕ+θ)^2)
                        *(
                          B*r^2*sin(ϕ+θ)^2 - E*(
                                                2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3^2*sin(θ)^2
                                               )
                          + C
                         )
                       )
           )
    return logarg
end # function D
