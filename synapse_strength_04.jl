# NOTES on version:
#   -01 - replaces C_pg and C_ph with fractions that are not caluclated in the
#           same way as other C parameters
#   -02 - adds in the option to use the Chavlis scheme or Santhakumar scheme for
#           defining the number of synpases
#   -03 - takes away the PP col and row, making the connectivity matrix only
#           used for the intenal connections.
#       - Adds C_scale to matrix
#       - Allows option of just using Chavlis' numbers, like a normal NMM
#   -04 - sets the MC-->GC connection to zero

using Statistics

function synapse_strength_04(C_scale, Model::String, mhc::String)
# - Calculates the synaptic strength, and hence the CONNECTIVITY PARAMETERS
#    of a Neural Mass Model.
# - Many of the strengths are estimated from Chavlis et al., 2017.
# - MC-HC connections are estimated from Larimer & Strowbridge, 2008.
# - This also seves as a folder for containing the paramaeters of Chavlis

###############################################################################
#################### Chavlis Charges ##########################################
###############################################################################

    #maximum conductances, Units = S
    #ampa
    g_A_max_pg = 0.8066 * 1e-9
    g_A_max_ph = 0.2400 * 1e-9
    g_A_max_gm = 0.5000 * 1e-9
    g_A_max_gb = 0.2100 * 1e-9
    g_A_max_mg = 0.1066 * 1e-9
    g_A_max_mb = 0.3500 * 1e-9
    #nmda
    g_N_max_pg = 0.8711 * 1e-9
    g_N_max_ph = 0.2760 * 1e-9
    g_N_max_gm = 0.5250 * 1e-9
    g_N_max_gb = 0.2310 * 1e-9
    g_N_max_mg = 0.1151 * 1e-9
    g_N_max_mb = 0.3850 * 1e-9
    #gaba
    g_G_max_bg = 14.000 * 1e-9
    g_G_max_hg = 0.1200 * 1e-9

    #conductance rise times, Units = s
    #ampa
    rise_A_pg = 0.10 * 1e-3
    rise_A_ph = 2.00 * 1e-3
    rise_A_gm = 0.50 * 1e-3
    rise_A_gb = 2.50 * 1e-3
    rise_A_mg = 0.10 * 1e-3
    rise_A_mb = 2.50 * 1e-3
    #nmda
    rise_N_pg = 0.33 * 1e-3
    rise_N_ph = 4.80 * 1e-3
    rise_N_gm = 4.00 * 1e-3
    rise_N_gb = 10.0 * 1e-3
    rise_N_mg = 0.33 * 1e-3
    rise_N_mb = 10.0 * 1e-3
    #gaba
    rise_G_bg = 0.90 * 1e-3
    rise_G_hg = 0.90 * 1e-3

    #conductance decay times, Units = s
    #ampa
    decay_A_pg = 2.500 * 1e-3
    decay_A_ph = 11.00 * 1e-3
    decay_A_gm = 6.210 * 1e-3
    decay_A_gb = 3.500 * 1e-3
    decay_A_mg = 2.500 * 1e-3
    decay_A_mb = 3.500 * 1e-3
    #nmda
    decay_N_pg = 50.00 * 1e-3
    decay_N_ph = 110.0 * 1e-3
    decay_N_gm = 100.0 * 1e-3
    decay_N_gb = 130.0 * 1e-3
    decay_N_mg = 50.00 * 1e-3
    decay_N_mb = 130.0 * 1e-3
    #gaba
    decay_G_bg = 6.800 * 1e-3
    decay_G_hg = 6.800 * 1e-3

    #resting potentials, Units = V
    v_rest_g = -85.0 * 1e-3
    v_rest_m = -64.0 * 1e-3
    v_rest_b = -52.0 * 1e-3
    v_rest_h = -59.0 * 1e-3
    #threshold potentials, Units = V
    v_th_g = -56.0 * 1e-3
    v_th_m = -42.0 * 1e-3
    v_th_b = -39.0 * 1e-3
    v_th_h = -50.0 * 1e-3
    #synapse reversal potential, Units = V
    E_ex = 0.0
    E_in = -70.0 * 1e-3
    #average membrane potential - important for the POST-synaptic neuron
    v_avg_g = mean([v_rest_g,v_th_g])
    v_avg_m = mean([v_rest_m,v_th_m])
    v_avg_b = mean([v_rest_b,v_th_b])
    v_avg_h = mean([v_rest_h,v_th_h])

    #arrays to make a loop easier
    g_max = [g_A_max_pg, g_A_max_ph, g_A_max_gm, g_A_max_gb, g_A_max_mg, g_A_max_mb, g_N_max_pg, g_N_max_ph, g_N_max_gm, g_N_max_gb, g_N_max_mg, g_N_max_mb, g_G_max_bg, g_G_max_hg]
    rise =  [ rise_A_pg,  rise_A_ph,  rise_A_gm,  rise_A_gb,  rise_A_mg,  rise_A_mb,  rise_N_pg,  rise_N_ph,  rise_N_gm,  rise_N_gb,  rise_N_mg,  rise_N_mb,  rise_G_bg,  rise_G_hg]
    decay = [decay_A_pg, decay_A_ph, decay_A_gm, decay_A_gb, decay_A_mg, decay_A_mb, decay_N_pg, decay_N_ph, decay_N_gm, decay_N_gb, decay_N_mg, decay_N_mb, decay_G_bg, decay_G_hg]
    v_avg = [   v_avg_g,    v_avg_h,    v_avg_m,    v_avg_b,    v_avg_g,    v_avg_b,    v_avg_g,    v_avg_h,    v_avg_m,    v_avg_b,    v_avg_g,    v_avg_b,    v_avg_g,    v_avg_g]
    E_syn = [E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_ex E_in E_in]

    #loop to calculate the charge
    q = zeros(1,length(g_max))
    for i=1:length(g_max)
        tau_tild = decay[i]*rise[i]/(decay[i]-rise[i])
        max_amp = (rise[i]/decay[i])^(tau_tild/decay[i])  -  (rise[i]/decay[i])^(tau_tild/rise[i])
        q[i] = (g_max[i]/max_amp) * (decay[i]-rise[i]) * (v_avg[i]-E_syn[i])

    end

    #synaptic receptor charges
    q_ampa_pg = q[1]
    q_ampa_ph = q[2]
    q_ampa_gm = q[3]
    q_ampa_gb = q[4]
    q_ampa_mg = q[5]
    q_ampa_mb = q[6]
    q_nmda_pg = q[7]
    q_nmda_ph = q[8]
    q_nmda_gm = q[9]
    q_nmda_gb = q[10]
    q_nmda_mg = q[11]
    q_nmda_mb = q[12]
    q_gaba_bg = q[13]
    q_gaba_hg = q[14]

    #in a model of membrane potential, NMDA and AMPA will not be modeled
    #individually, so average receptor charges to get an overall synapse charge
    qq = zeros(1,8)
    for n = 1:6
        # really want to know the rations between different receptor types here
        #qq[n] = q[n] + q[n+6]
        qq[n] = 0.5 * (q[n] + q[n+6])
    end
    qq[7] = q[13]
    qq[8] = q[14]

    #synaptic charges
    q_pg = qq[1]
    q_ph = qq[2]
    q_gm = qq[3]
    q_gb = qq[4]
    q_mg = qq[5]
    q_mb = qq[6]
    q_bg = qq[7]
    q_hg = qq[8]


###############################################################################
#################### Larimer and Strowbridge Charges ##########################
###############################################################################

# gl = inverse of the input resistance, Units = S
# ma = invese of the Maximum Amplitude of the double exponential formula
#      made up on the membrane constant (Decay) an the current decay
#      constant (Rise)
# A  = amplitude of the PSP, Units = V
# τm = membrane time constant (of post-syn. cell in connection)
# τd = current decay constant = PSP rise constant, Units = s

    τm_m = 22.0 * 1e-3 ########### where from
    τm_h = 12.0 * 1e-3 ########### where from
    #MC ---> HC
    gl_h  = 1/113e6    ########### where from
    A_mh  = 0.89 * 1e-3
    τd_mh = 3.06 * 1e-3
    T     = τm_h*τd_mh/(τm_h-τd_mh)
    ma_mh = 1/((τd_mh/τm_h)^(T/τm_h) - (τd_mh/τm_h)^(T/τd_mh))
    q_mh  = gl_h*A_mh*ma_mh*(τm_h -τd_mh)

    #HC ---> MC
    gl_m  = 1/120e6    ########### where from
    A_hm  = -0.42 * 1e-3
    τd_hm = 4.93 * 1e-3
    T     = τm_m*τd_hm/(τm_m-τd_hm)
    ma_hm = 1/((τd_hm/τm_m)^(T/τm_m) - (τd_hm/τm_m)^(T/τd_hm))
    q_hm  = gl_m*A_hm*ma_hm*(τm_m - τd_hm)








###############################################################################
######################## Matrix of Charges ####################################
###############################################################################

# - A matrix seems like the most intuitive way to order these conetions
# - There are 5 components to the model, so we will have a 5x5 Matrix
#                  P   G   M   B   H
#               P  -   -   -   -   -
#               G  -   -  g->m -   -
#               M  -   -   -   -   -
#               B  -   -   -   -   -
#               H  -   -   -   -   -
    Q = zeros(4,4) # matrix of raw charges
    #from GCs
    Q[1,2] = q_gm;   Q[1,3] = q_gb;
    #from MCs
    Q[2,1] = q_mg;   Q[2,3] = q_mb;   Q[2,4] = q_mh;
    #from BCs
    Q[3,1] = q_bg;
    #from HCs
    Q[4,1] = q_hg;   Q[4,2] = q_hm;
    Q = abs.(Q)
    Qrel = Q ./ maximum(Q)

###############################################################################
######################## Matrix of Synapses ###################################
###############################################################################

    #------------- Santhakumar ------------------------------------------------
    # - Can estimate this by looking at the number of synapses that a pre-syn cell
    #    cell makes with a post-syn cell. This is the 'divergence' in Santhakumar
    N_gm = 1.0
    N_gb = 1.0
    N_mg = 0.0
    N_mb = 1.0   #* mc_reduce
    N_mh = 2.0
    N_bg = 100.0
    N_hg = 160.0
    N_hm = 4.0
    Nsan = zeros(4,4)
    #from GCs
    Nsan[1,2] = N_gm;   Nsan[1,3] = N_gb;
    #from MCs
    Nsan[2,1] = N_mg;   Nsan[2,3] = N_mb;   Nsan[2,4] = N_mh;
    #from BCs
    Nsan[3,1] = N_bg;
    #from HCs
    Nsan[4,1] = N_hg;   Nsan[4,2] = N_hm;


    #--------- Chavlis -------------------------------------------------------
    Nchv = zeros(4,4)
    #from GC
    Nchv[1,2] = 16;    Nchv[1,3] = 20;
    #from MCs
    Nchv[2,1] = 0.0;   Nchv[2,3] = 100;
    #from BCs
    Nchv[3,1] = 20;
    #from HCs
    Nchv[4,1] = 16;
###############################################################################
######################## Cell proprtions ######################################
###############################################################################

    # These come from Chavlis et al. 2017
    M = zeros(4,4)
    M_p = 400  /2620
    M_g = 2000 /2620
    M_m = 80   /2620
    M_b = 100  /2620
    M_h = 40   /2620
    M = [M_g*ones(1,4);M_m*ones(1,4);M_b*ones(1,4);M_h*ones(1,4);]

###############################################################################
######################## Chavlis Connections ##################################
###############################################################################

# A total synapse matrix (NM) calculated by using the connectivity of the
# Chavlis model.
    # Total synapses
    NM = zeros(4,4)
    #from GC
    NM[1,2] = 32;   NM[1,3] = 2;
    #from MCs
    NM[2,1] = 0;    NM[2,3] = 8;
    #from BCs
    NM[3,1] = 2;
    #from HCs
    NM[4,1] = 16;


###############################################################################
######################## Synaptic Strengths ###################################
###############################################################################

# Strenghts are defined as
#                  C_xy = N_xy * Q_xy * M_x,
# where C_xy is the strength of the synpase between pre-syn population x and
# post-syn popultaion y.

    # Santhakumar numbers
    Csan      = Q.*Nsan.*M
    Csan[2,1] = 0.0
    Csan      = 100 * Csan ./ maximum(Csan) .* C_scale

    # Chavlis numbers
    Cchv      = Q.*NM
    Cchv[2,1] = 0.0
    Cchv      = 100 * Cchv ./ maximum(Cchv) .*C_scale

    # absorb the PP-->... connection parameters into the input



    if Model == "San"
        if   mhc == "mhcON"
            return Csan
        else mhc == "mhcOFF"
            Csan[2,4] = Csan[4,2] = 0.0
            return Csan
        end

    elseif Model =="Chg"
        return Q

    elseif Model == "Chv"
        return Cchv

    elseif Model == "ChvN"
        return Nchv .* C_scale

    else Model == "ChvNQ"
        return Nchv .* Qrel .* C_scale
    end

end
