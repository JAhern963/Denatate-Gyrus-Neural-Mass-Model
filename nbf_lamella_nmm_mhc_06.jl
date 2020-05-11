function nbf_lamella_nmm_mhc_06(YgA,YgN,Ym,Yb,Yh,n_gm_ampa,n_gm_nmda,rp,rgA,rgN,rm,rb,rh,t1,f,width,amp,n,C,Cpp,a_fac,tau_fac,U,plt)
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
    x_peaks = zeros(1, size(ti)[2])
    for i=1:size(ti)[2]
        x_peaks[i] = maximum(vx[ti[i]-1/(2*f) .< t .< ti[i]+1/(2*f)])
    end
    return x_peaks
end
peak_amps = zeros(5,10)
mc_peaks = peak_amps[1,:] = peaks(vm, ti, f, t)
bc_peaks = peak_amps[2,:] = peaks(vb, ti, f, t)
hc_peaks = peak_amps[3,:] = peaks(vh, ti, f, t)
yo_peaks = peak_amps[4,:] = peaks(yo, ti, f, t)
la_peaks = peak_amps[5,:] = peaks(lam, ti, f, t)
mc_fac_index = mean(mc_peaks[end-2:end]) ./ mc_peaks[1] # as defined in braganza 2020
bc_fac_index = mean(bc_peaks[end-2:end]) ./ bc_peaks[1]
yo_fac_index = mean(yo_peaks[end-2:end]) ./ yo_peaks[1]




######################## Plots ################################################

#=
print("\n")
print("mc_fac_index = ", mc_fac_index, ",   ")
print("bc_fac_index = ", bc_fac_index, ",   ")
print("signal_fac_index = ", yo_fac_index, "\n") =#

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
α = 0.4
plt1 = plot(t, vgA, legend=false, ylabel="Response (V)", xlabel="Time (s)",
    alpha=α, lc = :orange, label="GCa")
    plot!(t, vgN, lc=:orange, label="GCn", alpha=α, ls=:dot)
    plot!(t,  vm, lc=:green,  label="MC",  alpha=1)
    #plot!(ti[1:5], exp_data, marker=:circle, color=:green)
    plot!(t,  vb, lc=:blue,   label="BC",  alpha=α)
    plot!(t,  vh, lc=:purple, label="HI",  alpha=α)



### Plot of the reccording electode and its population components
plt2 = plot(t, C[2,1]*vm, ls=:dot, lc=:green, label="MC", legend=:topleft)
    plot!(t,          vp, ls=:dot, lc=:black, label="PP")
    plot!(t,  -C[3,1]*vb, ls=:dot, lc=:blue,  label="BC")
    plot!(t,  -C[4,1]*vh, ls=:dot, lc=:red,   label="HC")
    plot!(t,          yo, lc=:black,          label="Rec.",
    title="vp - C_bg.Vb - C_hg.Vh", legend=false)

### plot both of the outputs together
plt3 = plot(plt1, plt2, layout=(1,2), size=(1100,400))

if plt == 1
    PLT = plt1
    display(plt1)

elseif plt == 2
    PLT = plt3
    display(plt3)

elseif plt == 3
    PLT = 2
    display(plt2)

else plt == 0
end


return vp, vgA, vgN, vm, vb, vh, t, lam, peak_amps

end
