# NOTES on version
#       MHC_00 stands for Mossy Hilar Connections. This is the first vesion of
#           a Chavlis insipred NMM that includes the MC<-->HC connections that are
#           described in Larimer and Strowbridge, 2008.
#           Changes from the code lamella_nmm_chv_01:
#            - renamed synaptic gains
#            - all connectivity stuff in the synapse_strength_01.jl file
#            - more plots added
#       -01 - Change the rate constants
#       -02 - Puts the PP input inside the Sigmoid non-liniarity. Also has
#               individual FR curves
#       -03 - Adds in the Granule cell NMDA current
#       -04 - Redesigns PP input so that BCs and MCs also accept have it
#       -05 - Includes MC mass facilitation
#       -06 - Facilitation changed to a Tsodyk-Markram-like equation

using DifferentialEquations, Plots
include("input\\high_freq_protocol.jl")
include("Sigmoid.jl")
include("synapse_strength_04.jl")

#Parameter definitions and units
### Neural mass equations
# C     = Connectivity parameters (UNIT-LESS)
# n_gm  = _ampa/_nmda are the ratios of AMPARs/NMDARs for the GC-MC synapse
# r_x   = Rate constant for population X (HERTZ)
# tdd   = Time delay for PSP to travel one dendritc compartment (SECONDS)
#           = length of compartment / dednritic conduction speed
# Yi    = Excitatory (m and g), slow and fast inhibitory synaptic gains (VOLTS)
#
### Pulsed input
# t1    = time of first pulse (SECONDS)
# f     = frquency of pulses (HERTZ)
# amp   = amplitude of pulses (VOLTS) ???
# width = width of pulses (SECONDS), Note: Kats protocol used a 0.1ms wide pulse
# n     = number of pulses
#
### Sigmoid response function
# r     = Slope of sigmoid at x=v0 (/VOLTS.SECOND)
# e0    = 1/2 maximum firing rate from sigmoid (HERTZ)
# v0    = efective activity threshold for population (VOLTS)
#
### Facilitation
# a     = maximum facilitation current possible
# U     = constant controling the response of facilitation to firing rate
# tau_f = timescale for the facilitation
#
### Output
# ...from reccording electrode placed near proximal dendrites of granule cells
#    in the MML

######################## Parameters ###########################################
#synaptic gain
gain_tune = 1 * 1e-3
YgA = 3.0 * gain_tune ## x5 to get the facilitation of MCs
YgN = 3.0 * gain_tune
Ym  = 3.0 * gain_tune
Yb  = 20.0 * gain_tune
Yh  = 20.0 * gain_tune
#AMPAR/NMDAR ratios
n_gm_ampa = 1.0
n_gm_nmda = 0.5
#time/rate constants
#=
re  = 1 / 2e-3
rp  = re
rgA = re        + 1 / 40e-3
rgN = (re        + 1 / 40e-3) * 0.2
rm  = re        + 1 / 30e-3
rb  = 1 / 1e-3  + 1 / 10e-3
rh  = 1 / 7e-3  + 1 / 21e-3 =#

rp  = 1000   # 1ms
rgA = 1000   # 1ms
rgN = 200    # 5ms
rm  = 1000   # 1ms
rb  = 1500   # 1ms
rh  = 200    # 5ms

#=
rp = 1 / 2e-3
rg = 1 / 40e-3
rm = 1 / 30e-3
rb = 1 / 10e-3
rh = 1 / 21e-3

rp = 1 / 2e-3
rg = 100
rm = 100
rb = 1000
rh = 50
print([rgA rgN rm rb rh]) =#
#facilitation
a_fac   = 6
tau_fac = 0.05
U       = 20.0
#perforant path input
t1    = 2.0
f     = 50
width = 0.1 * 1e-3
amp   = 500
n     = 10
######################## Connectivity Parameters ##############################
Cpp       = [1.0, 0.125, 0.5, 0.0]
C_scale   = 1.0
C = synapse_strength_04(C_scale, "San", "mhcON")
C[4,3] = 1.0
#C[2,3] = 100  # mc excitattion to the bcs
#C[3,1] = 100  # to make BCs inhibit GCs more

######################## Changes from Notebook investigation ##################
#=
# mossy basket connection increases from ~3.5 to 20
C[2,3] = 50
# facilitation adjusted
a_fac = 4
tau_fac = 0.2
=#

######################## ODE solution #########################################
# code from this point onwards will be made into a function for easier use in
# the Juypter notebook enviroment


#system input
I = t-> high_freq_protocol(t,t1,f,amp,width,n)
protocol_start = t1 - (1/f)
protocol_end   = t1 + n*width + (n)*(1/f)
#pulse times
ti = t1.*ones(1,n)
for i = 2:n
    ti[i] = t1 + (i-1)*(1/f)
end

#sigmoid response function
e0 = 2.5
r  = 0.56 * 1/1e-3
v0 = 6.0  * 1e-3
Sg(x) = Sigmoid(x, e0, r, v0)
Sm(x) = Sigmoid(x, e0, r, v0)
Sb(x) = Sigmoid(x, e0, r, v0)          ### NOTE this needs to be thought about !
Sh(x) = Sigmoid(x, e0, r, v0)

p = [rp rgA rgN rm rb rh YgA YgN Ym Yb Yh tau_fac a_fac U I]
function chavlis_lamella!(du,u,p,t)
    vp, up, vgA, ugA, vgN, ugN, vm, um, vb, ub, vh, uh, lam = u
    rp, rgA, rgN, rm, rb, rh, YgA, YgN, Ym, Yb, Yh, tau_fac, a_fac, U, I = p
    #perforant path
    du[1] = up
    du[2] = 1.0*rp*I(t)  -  2*rp*up  -  rp^2*vp

    #granule AMPA
    gc_firing_rate = Sg(vp -(C[3,1]*vb + C[4,1]*vh))
    du[3] = ugA
    du[4] = YgA*rgA*gc_firing_rate - 2*rgA*ugA  - rgA^2*vgA
    #granule NMDA
    du[5] = ugN
    du[6] = YgN*rgN*gc_firing_rate - 2*rgN*ugN  - rgN^2*vgN

    #mossy
    mc_firing_rate = Sm(Cpp[2]*vp + C[1,2]*(n_gm_ampa*vgA + n_gm_nmda*vgN) - C[4,2]*vh)
    du[7] = um
    du[8] = Ym*rm*lam*mc_firing_rate - 2*rm*um - rm^2*vm
    mc_firing_rate = Sm(C[1,2]*n_gm_ampa*vgA  + Cpp[2]*vp)
    # line above means that only the GC AMPA and PP currents drive the facilitation

    #basket
    bc_firing_rate = Sb(Cpp[3]*vp + C[2,3]*vm + C[1,3]*vgA - C[4,3]*vh)
    du[9] = ub
    du[10] = Yb*rb*bc_firing_rate - 2*rb*ub - rb^2*vb

    #hilus interneuron
    hc_firing_rate = Sh(Cpp[4]*vp + C[2,4]*vm)
    du[11] = uh
    du[12] = Yh*rh*hc_firing_rate - 2*rh*uh - rh^2*vh

    #mossy mass facilitation
    du[13] = -lam/tau_fac + U*(a_fac - lam)*mc_firing_rate
end
#solve system
u0    = 1e-3*rand(1,13); u0[13]=0;
tspan = (0.0, protocol_end)
prob  = ODEProblem(chavlis_lamella!, u0, tspan, p)
sol   = solve(prob, Tsit5(), adaptive=false, dt=0.1e-3)

# define outputs
vp  = sol[1,:]
vgA = sol[3,:]
vgN = sol[5,:]
vm  = sol[7,:]
vb  = sol[9,:]
vh  = sol[11,:]
lam = sol[13,:]
t   = sol.t
yo  =  vp + C[2,1]*vm - (C[3,1]*vb + C[4,1]*vh)
#output =  v_m - (v_b + v_h)


# analise the output
function peaks(vx, ti, f, t)
    peaks = zeros(1, size(ti)[2])
    for i=1:size(ti)[2]
        peaks[i] = maximum(vx[ti[i]-1/(2*f) .< t .< ti[i]+1/(2*f)])
    end
    return peaks
end
peak_amps = zeros(4,10)
mc_peaks = peak_amps[1,:] = peaks(vm, ti, f, t)
bc_peaks = peak_amps[2,:] = peaks(vb, ti, f, t)
hc_peaks = peak_amps[3,:] = peaks(vh, ti, f, t)
yo_peaks = peak_amps[4,:] = peaks(yo, ti, f, t)
mc_fac_index = mean(mc_peaks[end-2:end]) ./ mc_peaks[1] # as defined in braganza 2020
bc_fac_index = mean(bc_peaks[end-2:end]) ./ bc_peaks[1]
yo_fac_index = mean(yo_peaks[end-2:end]) ./ yo_peaks[1]





######################## Plots ################################################

print("\n")
print("mc_fac_index = ", mc_fac_index, "\n")
print("bc_fac_index = ", bc_fac_index, "\n")
print("signal_fac_index = ", yo_fac_index, "\n")

# we only really care for the resposes after a stimulus input
disp_time = findall(x->(x>protocol_start)&(x<protocol_end), sol.t )
vp  = vp[disp_time]
vgA = vgA[disp_time]
vgN = vgN[disp_time]
vm  = vm[disp_time]
vb  = vb[disp_time]
vh  = vh[disp_time]
yo  = yo[disp_time]
lam = lam[disp_time]
t   = t[disp_time]

# experimental data (-ish) from Lysetskiy 2005
amp0 = maximum(vm[0.05.< t.< t1+1/(2*f)])
exp_data = [1.00,  1.92, 3.16, 3.78, 4.90] * amp0

### Plots of the average potentials of the populations
gr()
α = 0.2
plt1 = plot(t, vgA, legend=false, ylabel="Response (V)", xlabel="Time (s)",
    alpha=α, lc = :orange, label="GCa")
    plot!(t, vgN, lc=:orange, label="GCn", alpha=α, ls=:dot)
    plot!(t,  vm, lc=:green,  label="MC",  alpha=1)
    #plot!(ti[1:5], exp_data, marker=:circle, color=:green)
    plot!(t,  vb, lc=:blue,   label="BC",  alpha=α)
    plot!(t,  vh, lc=:purple, label="HI",  alpha=α)
display(plt1)

#=
### Plot of the reccording electode and its population components
plt2 = plot(t, C[2,1]*vm, ls=:dot, lc=:green, label="MC", legend=:topleft)
    plot!(t,          vp, ls=:dot, lc=:black, label="PP")
    plot!(t,  -C[3,1]*vb, ls=:dot, lc=:blue,  label="BC")
    plot!(t,  -C[4,1]*vh, ls=:dot, lc=:red,   label="HC")
    plot!(t,          yo, lc=:black,          label="Rec.",
    title="vp - C_bg.Vb - C_hg.Vh", legend=false)
display(plt2)
=#
