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

export symbintegral

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
    @register_symbolic ϵ1II(ξ, par1II)
    @register_symbolic ϵ1⊥(ξ, par1⊥)
    @register_symbolic ϵ2II(ξ, par2II)
    @register_symbolic ϵ2⊥(ξ, par2⊥)
    @register_symbolic ϵ3(ξ)

    ϵ1II(ξ, par1II) = ϵ2osc(ξ, par1II)
    ϵ1⊥(ξ, par1⊥) = ϵ2osc(ξ, par1⊥)
    ϵ2II(ξ, par2II) = ϵ2osc(ξ, par2II)
    ϵ2⊥(ξ, par2⊥) = ϵ2osc(ξ, par2⊥)
    ϵ3(ξ) = ϵ0(ξ)

    # Store trig values for ϕ
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # Wave numbers in the z direction

    @register_symbolic ρ1(r, ξ, c, par1⊥)
    @register_symbolic ρ2(r, ξ, c, par2⊥)
    @register_symbolic ρ3(r, ξ, c)
    @register_symbolic ρ̃1(r, ξ, c, par1II, par1⊥)
    @register_symbolic ρ̃2(r, ξ, c, par2II, par2⊥)

    # Derivatives of ρ would need to be defined separately, if later needed

    ρ1(r, ξ, c, par1⊥) = (r^2 + c^-2*ξ^2*ϵ1⊥(ξ, par1⊥))^(1/2)
    ρ2(r, ξ, c, par2⊥) = (r^2 + c^-2*ξ^2*ϵ2⊥(ξ, par2⊥))^(1/2)
    ρ3(r, ξ, c) = (r^2 + c^-2*ξ^2*ϵ3(ξ))^(1/2)
    ρ̃1(r, ξ, c, par1II, par1⊥) = (r^2 + (ϵ1II(ξ,par1II)/ϵ1⊥(ξ, par1⊥) - 1)*r^2 * cosϕ^2 + c^-2*ξ^2*ϵ1II(ξ,par1II))^(1/2)
    ρ̃2(r, ξ, c, par2II, par2⊥) = (r^2 + (ϵ2II(ξ,par2II)/ϵ2⊥(ξ, par2⊥) - 1)*(r*cosϕ*cos(θ) - r*sinϕ*sin(θ))^2 + c^-2*ξ^2*ϵ2II(ξ,par2II))^(1/2)

    # Gruesome algebraic expressions
    γ = (
         (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) * (
                             (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)+ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c)) - ϵ1⊥(ξ, par1⊥)*(ρ̃1(r,ξ,c,par1II,par1⊥)-ρ1(r,ξ,c,par1⊥))*(r^2*sinϕ^2-ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c))
                             /(ρ1(r,ξ,c,par1⊥)^2-r^2*sinϕ^2)
                            )
        *(
          (ϵ3(ξ)*ρ2(r,ξ,c,par2⊥)+ϵ2⊥(ξ, par2⊥)*ρ3(r,ξ,c)) - ϵ2⊥(ξ, par2⊥)*(ρ̃2(r,ξ,c,par2II,par2⊥)-ρ2(r,ξ,c,par2⊥)) * ((r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2-ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c))
          /(ρ2(r,ξ,c,par2⊥)^2 - (r*cosϕ*sin(θ)+r*sinϕ*cos(θ))^2)
         )
       )

    A = (
         (
         (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c))*(ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) - (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c))*(ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c))*exp(-2ρ3(r,ξ,c)*d)
         )
        *(
          (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)+ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c))*(ϵ3(ξ)*ρ2(r,ξ,c,par2⊥)+ϵ2⊥(ξ, par2⊥)*ρ3(r,ξ,c)) - (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)-ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c))*(ϵ3(ξ)*ρ2(r,ξ,c,par2⊥)-ϵ2⊥(ξ, par2⊥)*ρ3(r,ξ,c))*exp(-2ρ3(r,ξ,c)*d)
        )
        - (ρ̃1(r,ξ,c,par1II,par1⊥)-ρ1(r,ξ,c,par1⊥))*ϵ1⊥(ξ, par1⊥)/(ρ1(r,ξ,c,par1⊥)^2-r^2*sinϕ^2)
        *(
          (r^2*sinϕ^2 - ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c))*(ϵ3(ξ)*ρ2(r,ξ,c,par2⊥)+ϵ2⊥(ξ, par2⊥)*ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c))
          + 2*(ϵ2⊥(ξ, par2⊥)-ϵ3(ξ))*(
                          r^2*sinϕ^2 * (r^2*ρ1(r,ξ,c,par1⊥)-ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)^2) 
                          + ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)^2*(r^2 - 2*r^2*sinϕ^2 + ρ1(r,ξ,c,par1⊥)*ρ2(r,ξ,c,par2⊥))
                         )
          * exp(-2ρ3(r,ξ,c)*d) + (r^2*sinϕ^2+ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c))*(ϵ3(ξ)*ρ2(r,ξ,c,par2⊥)-ϵ2⊥(ξ, par2⊥)*ρ3(r,ξ,c))
          * (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c)) * exp(-4ρ3(r,ξ,c)*d)
        )
       )
    
    B = (
         (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)+ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) + 2*(ϵ1⊥(ξ, par1⊥)-ϵ3(ξ)) * (r^2*ρ2(r,ξ,c,par2⊥)-ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)^2-2*ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)^2)
         * exp(-2*ρ3(r,ξ,c)*d) + (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)-ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c)) * exp(-4*ρ3(r,ξ,c)*d)
        )
        + (ρ̃1(r,ξ,c,par1II,par1⊥)+ρ1(r,ξ,c,par1⊥))*ϵ1⊥(ξ, par1⊥)/(ρ1(r,ξ,c,par1⊥)^2-r^2*sinϕ^2)
        *(
          -1*(r^2*sinϕ^2-ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c))*(ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) + 2*(
                                                         r^2*sinϕ^2 *(ρ1(r,ξ,c,par1⊥)*ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)^2)
                                                         - ρ1(r,ξ,c,par1⊥)^2*ρ3(r,ξ,c)^2 + ρ1(r,ξ,c,par1⊥)*ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)^2
                                                        ) * exp(-2*ρ3(r,ξ,c)*d)
          - (r^2*sinϕ^2+ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c)) * exp(-4*ρ3(r,ξ,c)*d)
         )

    C = (ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)*(
               -1*(ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)+ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c)) + 2*ρ3(r,ξ,c)*(ϵ1⊥(ξ, par1⊥)-ϵ3(ξ))
               * (r^2+ρ1(r,ξ,c,par1⊥)*ρ2(r,ξ,c,par2⊥)) * exp(-2*ρ3(r,ξ,c)*d) + (ϵ3(ξ)*ρ1(r,ξ,c,par1⊥)-ϵ1⊥(ξ, par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c))
               * (ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c)) * exp(-4*ρ3(r,ξ,c)*d)
              )
    + (ρ̃1(r,ξ,c,par1II,par1⊥)-ρ1(r,ξ,c,par1⊥))*ϵ1⊥(ξ, par1⊥)/(ρ1(r,ξ,c,par1⊥)^2-r^2*sinϕ^2)*ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)
    *(
      (r^2*sinϕ^2-ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)+ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)+ρ3(r,ξ,c))
      + 2*ρ3(r,ξ,c)*(
              ρ1(r,ξ,c,par1⊥)^2*ρ2(r,ξ,c,par2⊥) + ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)^2 + r^2*sinϕ^2 * (ρ1(r,ξ,c,par1⊥)-ρ2(r,ξ,c,par2⊥))
             ) * exp(-2*ρ3(r,ξ,c)*d)
      - (r^2*sinϕ^2+ρ1(r,ξ,c,par1⊥)*ρ3(r,ξ,c)) * (ρ1(r,ξ,c,par1⊥)-ρ3(r,ξ,c)) * (ρ2(r,ξ,c,par2⊥)-ρ3(r,ξ,c)) * exp(-4*ρ3(r,ξ,c)*d)
     )
   )

    E = 4*ρ1(r,ξ,c,par1⊥)*ρ2(r,ξ,c,par2⊥)*ρ3(r,ξ,c)^2 * (ρ̃1(r,ξ,c,par1II,par1⊥)-ρ1(r,ξ,c,par1⊥))*ϵ1⊥(ξ, par1⊥)/(ρ1(r,ξ,c,par1⊥)^2-r^2*sinϕ^2) * exp(-2*ρ3(r,ξ,c)*d)

    D = (γ^-1 * (A - (ρ̃2(r,ξ,c,par2II,par2⊥)-ρ2(r,ξ,c,par2⊥))*ϵ2⊥(ξ, par2⊥)/(ρ2(r,ξ,c,par2⊥)^2-r^2*sin(ϕ+θ)^2)
                        *(
                          B*r^2*sin(ϕ+θ)^2 - E*(
                                                2*r^2*sinϕ*cos(θ)*sin(ϕ+θ) + ρ3(r,ξ,c)^2*sin(θ)^2
                                               )
                          + C
                         )
                       )
           )
    rlnD = r*log(D)
    lnD = log(D)

    # Define derivative of lnD with respect to θ
    ∂θ = Differential(θ)
    r∂lnD∂θ = expand_derivatives(∂θ(rlnD))
    ∂lnD∂θ = expand_derivatives(∂θ(lnD))

    # Build Julia function for ∂lnD∂θ
#    function ∂lnD∂θfn(θval, rval, ϕval, ξval, par1IIval, par1⊥val, par2IIval, par2⊥val, dval, cval)
#        return build_function(∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
#                              expression = Val{false})
#    end # function ∂lnD∂θfn
    r∂lnD∂θexpr = build_function(r∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                                   expression = Val{true})
    r∂lnD∂θfn = build_function(r∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                                   expression = Val{false})
    ∂lnD∂θexpr = build_function(∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                                expression = Val{true})
    ∂lnD∂θfn = build_function(∂lnD∂θ, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                                expression = Val{false})


    #println(expand_derivatives(∂lnD∂θ))

#    ## Try evaluating ∂lnD∂θ at some point
#    valuedict = Dict([θ => θval, r => 50119u"m^-1", ϕ => 0.58, ξ => 0u"Hz2π", par1II => calciteII, 
#                      par1⊥ => calcite⊥, par2II => BaTiO3II, par2⊥ => BaTiO3⊥, 
#                      d => dval, c => c_0])
#
#    println("∂lnD∂θ = $(substitute(∂lnD∂θ, valuedict))")

"Attempts to find symbolic integrals over ϕ and r"
function symbintegral()
    intega = integrate(A, ϕ, symbolic=true, detailed=false)
    #integϕ = integrate(r∂lnD∂θexpr, ϕ, symbolic=true, detailed=false)
    #integr = integrate(integϕ, r, symbolic=true, detailed=false)

    integrfcn = build_function(intega, θ, r, ϕ, ξ, par1II, par1⊥, par2II, par2⊥, d, c,
                              expression = Val{false})
    return integrfcn
end # function symbintegral

"Returns the torque per unit area of the two-plate system with dielectric functions 
ϵ1II(ξ), ϵ1⊥(ξ), ϵ2II(ξ) and ϵ2⊥(ξ), ϵ3(ξ), separation dval, principal axes angle θval and temperature T. Uses a trapeziodal rule if trapz=true. A value rcutoff = 0 m^-1 will take rcutoff=d^-1"
function twoplateM(dval, θval, T; ξcutoff = 1e18u"Hz2π", rcutoff = 0u"m^-1", rtol = 1e-2, atolunitful=1e-19u"J*m^-2", integmethod="trapz", trapzrdiv=100, trapzϕdiv=100, ξcutofffactor=1e-5)
    # r is the transverse wavenumber

    if !(typeof(T) <: Unitful.Temperature)
        T = T*u"K"
        println("Temperature was not given in temperature Unitful type, assuming Kelvin")
        end # if

    atol = ustrip(u"m^-2", atolunitful/(k_B*T/(4*π^2)))

    # Matsubara sum #
    M = 0u"J*m^-2"
    Merror = 0u"J*m^-2"
    ξn = 0u"Hz2π"
    #ξn = uconvert(u"Hz2π", k_B*T/ħ)


    if rcutoff == 0u"m^-1"
        rcutoff = dval^-1
    end # if
    rcutoffnounit = ustrip(u"m^-1", rcutoff) # Strip units for hcubature,
    #has to match inserted length unit below
    
    ξloopcounter = 0 # Counts subsequenct frequencies that haven't added to the result

    while ξn < ξcutoff
#        integ = rϕ -> rϕ[1] * substitute(∂lnD∂θ, Dict([θ => θval, r => rϕ[1]*u"m^-1", ϕ => rϕ[2], 
#                                                    ξ => ξn, par1II => calciteII, 
#                                                    par1⊥ => calcite⊥, par2II => BaTiO3II,
#                                                    par2⊥ => BaTiO3⊥, d => dval, c => c_0]))
        function integtrapz(r, ϕ)
            #println("Calculating torque for r = $(rϕ[1]), ξ=$ξn, and ϕ = $(rϕ[2])")
            integrandval = ∂lnD∂θfn(θval, r, ϕ, ξn, calciteII, calcite⊥,
                                             BaTiO3II, BaTiO3⊥, dval, c_0)
            yield() # Lousy attempt to allow for interruption
            if isnan(integrandval)
                println("NaN at r = $(r), ξ=$ξn, and ϕ = $(ϕ)")
                throw(ErrorException("NaN at r = $(r), ξ=$ξn, and ϕ = $(ϕ)"))
            end # if
            return integrandval
        end # function integ

        function integ(rϕ)
            #println("Calculating torque for r = $(rϕ[1]), ξ=$ξn, and ϕ = $(rϕ[2])")
            integrandval = r∂lnD∂θfn(θval, rϕ[1]*u"m^-1", rϕ[2], ξn, calciteII, calcite⊥,
                                             BaTiO3II, BaTiO3⊥, dval, c_0)
            yield() # Lousy attempt to allow for interruption
            if isnan(integrandval)
                println("NaN at r = $(rϕ[1]), ξ=$ξn, and ϕ = $(rϕ[2])")
                throw(ErrorException("NaN at r = $(rϕ[1]), ξ=$ξn, and ϕ = $(rϕ[2])"))
            end # if
            return ustrip(u"m^-1", integrandval)
        end # function integ

        # Split integration into two domains in ϕ, to avoid divisions by zero
        # at sin(ϕ) = ± 1
        if integmethod == "trapz"
            η = 1e-5u"m^-1"
            rr = range(η, rcutoff, length=trapzrdiv)
            ϕϕ = range(0, 2π, length=trapzϕdiv)
            mm = [integtrapz(r, ϕ) for r in rr, ϕ in ϕϕ] # Lacks prefactor, also r dr missing
            phiint = uconvert(u"aJ", k_B*T/(4*π^2)) * trapz(ϕϕ, mm, Val(2))
            phiint[phiint .< 1e-5u"aJ"] .= 0u"aJ"
            phiint = phiint .* rr # Add factor r
            println("Integrated over ϕ for ξn = $ξn")
            #display(stack((reverse(phiint), reverse(rr))))
            rint = trapz(rr, phiint)
            
            if rint < ξcutofffactor*M
                ξloopcounter += 1
            end # if
            M += rint
        elseif integmethod == "mixed"
            η = 1e-5u"m^-1"
            rr = range(η, rcutoff, length=trapzrdiv)

            integmix = [ϕ -> integtrapz(r, ϕ) for r in rr] # Lacks prefactor, also r dr missing

            atoleffective = ustrip.(u"m^2", atol./(rr*rcutoff))
            ϕintervals = (-π/2, π/4, π/2, 5*π/4, 3*π/2) 
            # Above probably depend on θ and so should be calculated
            phiint = []
            phiint .+= [uconvert(u"aJ", k_B*T/(4*π^2)) * quadgk(integmix[i], ϕintervals...,
                                                                rtol=rtol,
                                                                atol=atoleffective[i])[1]
                          for (i, r) in enumerate(rr)]
            phiint[phiint .< 1e-5u"aJ"] .= 0u"aJ"
            phiint = phiint .* rr # Add factor r
            println("Integrated over ϕ for ξn = $ξn")
            display(stack((reverse(phiint), reverse(rr))))
            rint = trapz(rr, phiint)
            
            if rint < ξcutofffactor*M
                ξloopcounter += 1
            end # if
            M += rint
        elseif integmethod == "hcubature"
            Marray = (k_B*T/(4*π^2) .* hcubature(integ, (0, -π/2), (rcutoffnounit, π/2), rtol=rtol)
                                                 #atol=atol)
                      .*u"m^-2")
            M += Marray[1]
            Merror += Marray[2]
            Marray = (k_B*T/(4*π^2) .* hcubature(integ, (0, π/2), (rcutoffnounit, 3π/2), rtol=rtol)
                                                 #atol=atol)
                      .*u"m^-2")
            M += Marray[1]
            Merror += Marray[2]
        else
            throw(ErrorException("Integration method $integmethod not recognized"))
        end # if
#        Marray = (k_B*T/(4*π^2) .* hcubature(integ, (0, π/2), (rcutoff, 3π/2), rtol=rtol)
#                                             #atol=atol)
#                  .*u"m^-2")
#        M += Marray[1]
#        Merror += Marray[2]
        # Unit of r dr is missing

        if ξn == 0 # Halve first term
            M = M/2
        end # if

        if ξloopcounter > 2
            println("Last few Matsubara terms has not contributed significantly to the torque, breaking loop")
            break
        end # if

        println("M contribution from ξn = $ξn calculated, current torque M = $M")
        ξn += uconvert(u"THz2π", k_B*T/ħ)
    end # while

    return M
end # function twoplateM

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

"Plots the torque contribution for twoplateM, over r and ϕ for a given ξ. If rfac is true,
the integration measure factor r is included."
function Mplot(T, ξ, par1II, par1⊥, par2II, par2⊥, dval, θval; rmin = -4, rmax=15, ϕmin=0, ϕmax=2π, ϕstep=0.1, rstep=0.5, rfac=true)
    rrange = (10. .^ (rmin:rstep:rmax))u"m^-1"
    ϕrange = ϕmin:ϕstep:ϕmax
    MM = zeros(length(rrange), length(ϕrange))
    if rfac
        MM = MM .* u"m^-1"
    end # if
    for i in 1:length(rrange)
        for j in 1:length(ϕrange)
            try
                if rfac
                    MM[i, j] = r∂lnD∂θfn(θval, rrange[i], ϕrange[j], ξ, par1II, par1⊥, par2II, par2⊥, dval, c_0)
                else
                    MM[i, j] = ∂lnD∂θfn(θval, rrange[i], ϕrange[j], ξ, par1II, par1⊥, par2II, par2⊥, dval, c_0)
                end # if
            catch e
                println("Error at r = $(rrange[i]) and ϕ = $(ϕrange[j])")
                rethrow()
            end # try
        end # for
    end # for

    MM = MM.* (k_B*T/(4*π^2)) # Now adds prefactor, such that the integral can be estimated 
                          # directly from the plot

    if rfac
        title = "k_B T/4π² r∂lnD∂θ(r, ϕ) for ξ = $ξ"
    else
        title = "k_B T/4π² ∂lnD∂θ(r, ϕ) for ξ = $ξ"
    end # if

    ## Set the unit to something reasonable, so that low valued don't mess with
    #the heatmap etc
    
    if rfac
        MM = uconvert.(u"aJ*m^-1", MM)
    else
        MM = uconvert.(u"aJ", MM)
    end # if
    Mmax = maximum(MM)
    Mmin = minimum(MM)
    println("Mmax = $Mmax, Mmin = $Mmin")
    #Mstep = (Mmax - Mmin)/5
    #colorbar_ticks=(collect(Mmax:Mstep:Mmin), Printf.format.(Ref(Printf.format"%.3E"), Mmax:Mstep:Mmin))
    #println(colorbar_ticks)
    
    l = @layout [a{0.98h} ; b]
    p = heatmap(ϕrange, rrange, MM, xlabel="ϕ", ylabel="r", title=title,
            yscale=:log10)
    p2 = plot(axis=([], false))
    ftr = text("θ=$θval, d=$dval, ϕstep=$ϕstep, rstep=$rstep", :black, :right, 8)
    annotate!(1, 0.5, ftr)
    plot(p, p2, layout=l)
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
