%%
% Default parameter values for Kumar et al., 2023
%
% Takes in the firing rate of the background, whether pre- (0) or post- (1)
% damage, and the number of populations (3: PN, PV, and SOM; 4: same + VIP)
%%
function [params] = A1_params(r_bg,damaged,Npop)

if Npop ~= 3 && Npop ~=4
   error('Code only tested for 3 and 4 neuronal populations'); 
end

%% Numerical parameters for spiking sims
params.T = 2000;  %Default: 5000; % total time
params.dt = 0.01;
params.NT = round(params.T/params.dt); % total number of time steps
params.recstart = 1000; % time at which to start analysis, to remove transient

%% Population parameters
% Save Npop in params structure
params.Npop = Npop;

% Set the number of neurons in each population
% Order: [e, pv, som, vip]
params.Ncells = zeros(params.Npop,1); 
params.Ncells(1)= 5000;
params.Ncells(2:Npop) = 520;
params.Ntot = sum(params.Ncells);

% start index of each population (ends with Ntot)
params.pinds = cumsum([0;params.Ncells])+1;   

% Preallocate rates and set the max number of spikes for the spiking sims
params.rates = zeros(params.Ntot,1);
params.maxns = round(params.T*params.Ntot*0.05);

%% Connectivity parameters

% connection probabilities
% Columns are pre-synaptic populations [PN, PV, SOM, VIP]
% Rows are post-synaptic populations [PN, PV, SOM, VIP]
params.p = [0.03 0.10 0.10 0.0;
    0.05 0.10 0.07 0.0;
    0.05 0.00 0.00 0.1;
    0.05 0.15 0.05 0.0];
% Extract the necessary submatric for Npop = 3
params.p = params.p(1:Npop,1:Npop);
% Note that we take the transpose of p for calculations in the code
params.p = params.p';

% Connection weights
params.w = 0.6;
params.g = 3;
params.W = params.w*[1, -params.g, -params.g, -params.g;
                     1, -params.g, -params.g, -params.g;
                     1, -params.g, -params.g, -params.g;
                     1, -params.g, -params.g, -params.g];
params.W = params.W(1:Npop,1:Npop);         
params.gSyn = params.W';   
params.J = params.gSyn;

%% Feedforward connections
% probability of an external connection
params.pEext = 0.1;
params.pPext = 0.08;
params.pSext = 0.0;
params.pVext = 0.08;

% Number of external connections (assume the number of external inputs  
% matches the number of excitatory cells)
params.N_ext = [params.pEext  params.pPext params.pSext params.pVext]*params.Ncells(1);
params.N_ext = params.N_ext(1:Npop);

% total strength of feedforward connections
params.J_ext = params.w*params.N_ext;
params.J_sigma_ext = params.w^2*params.N_ext;

%% Neuron-specific parameters (all cells have the same intrinsic properties 

% cell parameters, one value for each population
params.EL =  -65*ones(1,Npop); %leak voltage, mV
params.vth = -50*ones(1,Npop); %threshold voltage, mV
params.vre = -65*ones(1,Npop); %reset voltage, mV

params.tau_m = 10*ones(1,Npop); %Cm./gL, membrane time constant (ms);
params.tau_s = 0.5*ones(1,Npop); % synaptic time constant (ms)
params.tau_ref = 2*ones(1,Npop); %refractory period (ms)

%% Mean and standard deviation of baseline input (default and damaged)
params.mu_bg = params.J_ext.*params.tau_m*(r_bg*1e-3);
params.sigma_bg = sqrt(params.J_sigma_ext.*params.tau_m*(r_bg*1e-3));
params.sigma_fixed = sqrt(10*ones(1,Npop));

if damaged == 0 % i.e., not damaged
    params.bg_damage = [1 1 1 1];
    params.stim_damage = [1 1 1 1];
    params.recov = [0 0 0 0];
else % example params for the damaged condition
    
    params.bg_damage = [0.5 0.5 0 0.5];
    if Npop == 3
        params.stim_damage = [0.1 0.1 0];
        params.recov = [3 2 -5];
    else
        params.stim_damage = [0.1625 0.05 0 0.5];
        params.recov = [2.5 0 0 1.25];
    end
end
% Extract the relevant values
params.bg_damage = params.bg_damage(1:Npop);
params.stim_damage = params.stim_damage(1:Npop);
params.recov = params.recov(1:Npop);

%% Mean and standard deviation of stimulus
% The firing rate of the stimulus at different intensities
nu_stim = [0 2 4 8];
params.num_stims = length(nu_stim);

% The feedforward stimulus is excitatory, so the jump in mu is given by
% the following
for ii = 1:params.num_stims
    params.mu_stim(ii,:) = params.J_ext.*params.tau_m*(nu_stim(ii)*1e-3);
    params.sigma_stim(ii,:) = sqrt(params.J_sigma_ext.*params.tau_m*(nu_stim(ii)*1e-3));
end

%% Corresponding parameters for theory (diffusion approx, units: seconds)
params.alpha = sqrt(2)*abs(zeta(0.5));
params.V_th = 15;
params.V_r = 0;

% in-degrees    
params.I = params.p'.*params.Ncells';
params.J_theory = params.I.*params.W;
params.J_sigma_theory = (params.I.*params.W.*params.W);

params.tau_m_theory = params.tau_m(1)*1e-3; %Cm./gL, membrane time constant (ms);
params.tau_s_theory =  params.tau_s(1)*1e-3; % synaptic time constant (ms)
params.tau_ref_theory = params.tau_ref(1)*1e-3; %refractory period (ms)


end

