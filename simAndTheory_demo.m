%%
% This code provides a demo for the spiking simulation and mean field 
% theory for Kumar et al., 2023
%
% The spiking simulation uses mex code and a par-for loop
% The simulations run for 2 seconds (of simulation time), which is shorter 
% than the default 5 second, but finishes faster (~30 seconds depending on
% computer)
%%
clear; clc; close all;

restoredefaultpath;
folder = fileparts(which('simAndTheory_demo.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Load the parameters (for both the spiking sims and the theory)
r_bg = 3; % firing rates of background input
damaged = 0; % pre-damaged: 0; post-damaged: 1
Npop = 4; % number of populations
params = A1_params(r_bg,damaged,Npop);

%% Run the spiking simulations

% Compile the mex code (do after each edit to the mex file)
% Run the following if running mex for the first time: run mex -setup -c
mex Mex_Functions/LIF_diff_approx_mex.c

% seed the random number generator
% rng('shuffle'); % use this for new random numbers
rng(314); % use this for consistency/reproducibility 

% Get the currently used random seed
s = rng(); 
random_seed = s.Seed;

% Preallocate
rates_trial = zeros(params.num_stims,sum(params.Ncells));
times = zeros(params.num_stims,params.maxns);
tinds = zeros(params.num_stims,params.maxns);

% Loop through the stimulus strengths (none, low, med, and high)
tic;
parfor par_loop = 1:params.num_stims
    % Define the feedforward input (mean and variance) accounting for possible damage
    mu_ext = params.mu_stim(par_loop,:).*params.stim_damage + ...
        params.mu_bg.*params.bg_damage+params.recov;
    var_ext = sqrt(params.sigma_stim(par_loop,:).^2.*params.stim_damage+...
        params.sigma_bg.^2.*params.bg_damage+params.sigma_fixed.^2);

    % Create the connectivity matrix
    rng(random_seed); % Comment out for different connectivity matrices across stims
    [wind,wipost,wstr] = ...
        Diff_approx_gen_weights(params.Ncells,params.p,params.J);
    
    % offset the index by one for the MEX code (C starts indexing at 0)
    wind_mex = int32(wind-1);
    wipost_mex = int32(wipost-1);
    pinds_mex = int32(params.pinds-1);
         
    % Run the simulation
    [~,times(par_loop,:),tinds(par_loop,:)] = ...
        LIF_diff_approx_mex(params.T, params.NT, params.recstart,...
        params.Npop, params.Ntot, params.Ncells, params.dt, params.rates, params.p, params.J,...
        params.tau_s, params.tau_m, params.EL, params.vth, params.vre, params.tau_ref,...
        wind_mex, wipost_mex, wstr, pinds_mex,random_seed,mu_ext,var_ext);   
end
toc;

% Process the spiking data (i.e., get average firing rates)
spikeRatesAve = zeros(params.Npop,params.num_stims);
for ii = 1:params.num_stims
    spikeRatesAve(:,ii) = spikeAnalysisFn(tinds(ii,:),times(ii,:),params);
end
p = polyfit([0 1 2 3],spikeRatesAve(1,:),1);
fprintf('Spiking: The gain estimate for E is %.2f \n',p(1))

%% Run the mean field theory
[firing_rates_sol,nan_warning,warning_notice, eig_values, d_min, real_eig_max]...
    = stim_loop_fn(params,params.bg_damage,params.stim_damage,params.recov);
p = polyfit([0 1 2 3],firing_rates_sol(1,:),1);
fprintf('Mean field: The gain estimate for E is %.2f \n',p(1))

%% Create a plot comparing the firing rates
color_scheme =[59, 57, 60; 164, 71, 105;181, 117, 51; 107, 159, 165]/255;

figure(1); clf; hold on;
hLeg = [];
for ii = 1:params.Npop
    hLeg(ii) = plot(spikeRatesAve(ii,:),'.-','markersize',15,'color',color_scheme(ii,:));
    plot(firing_rates_sol(ii,:),'*','markersize',16,'color',color_scheme(ii,:))
end

set(gca,'fontsize',16)
popNames = {'PN','PV','SOM','VIP'};
legend(hLeg,popNames(1:params.Npop))
xlabel('Stimulus Strength')
ylabel('Firing Rate (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})
box off

%% Create raster plots
figure(2); clf; hold on;
titles = {'None','Low','Med','High'};
for ii = 1:params.num_stims
    subplot(1,params.num_stims,ii)
    plot_raster(times(ii,:), tinds(ii,:), params.Ntot, params.Npop, ...
        params.Ncells, params.pinds, 4)
    xlim([1 2])    
    if ii >= 2
        legend off
    else
        ylabel('Neuron Index')
    end
    title(titles{ii})
end
