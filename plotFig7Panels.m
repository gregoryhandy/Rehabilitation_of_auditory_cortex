%%
% Reproduces the modeling panels in Kumar et al., 2023
%
% Uses compressed spiking datasets for the 4-population spiking simulations
% (panel b and c) 
% Reruns the code for the mean field theory for panel c (i.e., simulations 
% vs theory firing rates)
% Uses the compressed datasets corresponding to the parameter search for 
% the mean field model (panels d, e, and f)
%
% Code is written in blocks but should be run sequentially
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('plotFig7Panels.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Plot Fig 7b (raster plots)
% Load the spiking data
spikingData = load('Sim_Data/fourPopDefault.mat');

% Create the raster plots
titles = {'None','Low','Med','High'};
figure(1); clf;
for jj = 1:4
    subplot(1,4,jj)
    plot_raster(spikingData.times(jj,:), spikingData.tinds(jj,:), spikingData.params.Ntot,...
        spikingData.params.Npop, spikingData.params.Ncells, spikingData.params.pinds, 4)
    xlim([3 4])
    if jj >= 2 
        legend off
    else
        ylabel('Neuron Index')
    end
    title(titles{jj})
end

%% Plot Fig 7c (simulations vs. mean field firing rates)

% Process the spiking data (i.e., get average firing rates)
spikeRatesAve = zeros(4,4);
for ii = 1:4
    spikeRatesAve(:,ii) = spikeAnalysisFn(spikingData.tinds(ii,:),...
        spikingData.times(ii,:),spikingData.params);
end
p = polyfit([0 1 2 3],spikeRatesAve(1,:),1);
fprintf('Spiking: The gain estimate for E is %.2f \n',p(1))

% Load parameters and run the mean field theory
params = A1_params(3,0,4);
[firing_rates_sol,nan_warning,warning_notice, eig_values, d_min, real_eig_max]...
    = stim_loop_fn(params,params.bg_damage,params.stim_damage,params.recov);

% Calculate gain
p = polyfit([0 1 2 3],firing_rates_sol(1,:),1);
fprintf('Mean field: The gain estimate for E is %.2f \n',p(1))

% Create the plot 
figure(2); clf; hold on;
h = [];
color_scheme =[59, 57, 60; 164, 71, 105;181, 117, 51; 107, 159, 165]/255;
for jj = 1:params.Npop
    h(jj) = plot(1:4,spikeRatesAve(jj,:),'.-','markersize',15,'color',color_scheme(jj,:));  
    plot(1:4,firing_rates_sol(jj,:),'*','markersize',16,'color',color_scheme(jj,:))  
end
set(gca,'fontsize',16)
legend(h,{'PN','PV','SOM','VIP'})
xlabel('Stimulus Strength')
ylabel('Firing Rate (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})
box off

%% Plot Fig 7d (firing rate curves for the different parameters)

% Load the condensed data set for the parameter sweep
load('Sim_Data/paramSweepFourPop.mat')
num_samples = length(improved_gain_indices);

figure(3); clf; hold on;
for jj = 1:params.Npop
    stim_levels = [0 1 2 3];
    for ii = 1:num_samples
        plot(stim_levels,squeeze(firing_rates_sorted(ii,jj,:)),'-.',...
            'color',[color_scheme(jj,:), 1],'linewidth',0.5);
    end
end

h = [];
for jj = 1:params.Npop
    mean_est = mean(squeeze(firing_rates_sorted(1:num_samples,jj,:)));
    h(jj) = plot(stim_levels,mean_est,'color',color_scheme(jj,:),'linewidth',4);
    
    set(gca,'fontsize',16)
    xticks([0 1 2 3])
    xticklabels({'None','Low','Med','High'})
    xlabel('Stimulus Strength')
    ylabel('Firing rate (Hz)')
end

legend(h,{'PN','PV','SOM', 'VIP'})
box off

%% Plot Fig 7e (distribution of successful recovery currents)

plot_titles = {'PN','PV','SOM','VIP'};
popsToPlot = [1, 2, 4];

figure(4); clf; hold on;
% Set the edges for the histogram
histEdges = linspace(0,5,5);
histEdges = [histEdges - diff(histEdges(1:2))/2, histEdges(end)+diff(histEdges(1:2))/2];
for ii = 1:3
    subplot(1,3,ii) 
    histogram(recov_amounts_sorted(improved_gain_indices,popsToPlot(ii)), histEdges,...
        'Normalization','pdf','Facecolor',color_scheme(popsToPlot(ii),:),'Edgecolor',color_scheme(popsToPlot(ii),:))
    
    set(gca,'fontsize',16)
    xlabel('I_{recov}')
    title(plot_titles{popsToPlot(ii)})
    box off
end

%% Plot Fig 7f (box plots of the input into SOM cells)

% Get the input currents
params = A1_params(3,1,4);
mu_SOM = zeros(4,length(improved_gain_indices));
for ii = 1:length(improved_gain_indices)
    for jj = 1:4
        mu_ext = params.mu_stim(jj,:).*damage_amounts_sorted(ii,:) + ...
        params.mu_bg.*params.bg_damage+recov_amounts_sorted(ii,:);
        mu_all = params.J_theory*firing_rates_sorted(ii,:,jj)'*params.tau_m_theory+mu_ext';
        mu_SOM(jj,ii) = mu_all(3);
    end
end

% Create the plot
figure(5); clf; hold on;
boxplot(mu_SOM','Notch','on','Whisker',10,'colors',color_scheme(3,:))
default_input =  [4.1739, 7.8299, 9.1286, 10.9058];
for ii = 1:4
    h = plot(ii, default_input(ii),'k.','markersize',16);
    set(gca,'fontsize',16)
    ylim_temp = ylim;
    ylim([ylim_temp(1),ylim_temp(2)+1])
    ylabel('Ave. Synaptic Input to SOM (mV)')
    xlabel('Stimulus Intensity')
    xticks([1 2 3 4])
    xticklabels({'None', 'Low','Med','High'})
end

% Fill in the box plots
gg = findobj(gca,'Tag','Box');
boxplot_lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
for j=1:length(gg)
    patch(get(gg(j),'XData'),get(gg(j),'YData'),color_scheme(3,:),'FaceAlpha',0.5);
    set(boxplot_lines, 'Color', color_scheme(3,:));
end

legend(h,'Control (Undamaged)')
box off
