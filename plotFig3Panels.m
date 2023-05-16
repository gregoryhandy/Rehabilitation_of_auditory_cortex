%%
% Reproduces the modeling panels in Fig 3 for Kumar et al., 2023
%
% Uses compressed spiking datasets for the 3-population spiking simulations
% (panel b, c, and d) 
% Reruns the code for the mean field theory for panel c (i.e., simulations 
% vs theory firing rates)
% Uses the compressed datasets corresponding to the parameter search for 
% the mean field model (panels e, and f)
%
% Code is written in blocks but should be run sequentially
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('plotFig3Panels.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Plot 3b (raster plots)
% Load the spiking data
spikingData = load('Sim_Data/threePopDefault.mat');

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

%% Plot Fig 3c (simulations vs. mean field firing rates)

% Process the spiking data (i.e., get average firing rates)
spikeRatesAve = zeros(3,4);
for ii = 1:4
    spikeRatesAve(:,ii) = spikeAnalysisFn(spikingData.tinds(ii,:),...
        spikingData.times(ii,:),spikingData.params);
end
p = polyfit([0 1 2 3],spikeRatesAve(1,:),1);
fprintf('Spiking: The gain estimate for E is %.2f \n',p(1))

% Load parameters and run the mean field theory
params = A1_params(3,0,3);
[firing_rates_sol,nan_warning,warning_notice, eig_values, d_min, real_eig_max]...
    = stim_loop_fn(params,params.bg_damage,params.stim_damage,params.recov);

% Calculate gain
p = polyfit([0 1 2 3],firing_rates_sol(1,:),1);
fprintf('Mean field: The gain estimate for E is %.2f \n',p(1))

% Create the plot 
figure(2); clf; hold on;
h = [];
color_scheme =[59, 57, 60; 164, 71, 105;181, 117, 51]/255;
for jj = 1:params.Npop
    h(jj) = plot(1:4,spikeRatesAve(jj,:),'.-','markersize',15,'color',color_scheme(jj,:));  
    plot(1:4,firing_rates_sol(jj,:),'*','markersize',16,'color',color_scheme(jj,:))  
end
set(gca,'fontsize',16)
legend(h,{'PN','PV','SOM'})
xlabel('Stimulus Strength')
ylabel('Firing Rate (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})
box off

%% Plot Fig 1d (example raster plots for the insets) 

condsToPlot = [1 4];
figure(3); clf;
for outerLoop = 1:2
    if outerLoop == 1
        spikingDataEx = load('Sim_Data/threePopUnstableEx.mat');
    else
        spikingDataEx = load('Sim_Data/threePopStableEx.mat');
    end
    
    for jj = 1:2
        subplot(2,2,jj+2*(outerLoop-1))
        plot_raster(spikingDataEx.times(condsToPlot(jj),:), spikingDataEx.tinds(condsToPlot(jj),:), ...
            spikingDataEx.params.Ntot, spikingDataEx.params.Npop, spikingDataEx.params.Ncells,...
            spikingDataEx.params.pinds, 4)
        xlim([3 4])
        if jj >= 2 || outerLoop == 1
            legend off
        end
        
        if jj == 1
            ylabel('Neuron Index')
        end
        
        if outerLoop == 1
            title(titles{jj})
        end
    end
end

%% Plot 3d (firing rate curves for the different parameters)

% Load the condensed data set
paramSweepResults = load('Sim_Data/paramSweepThreePop.mat');
num_samples = length(paramSweepResults.improved_gain_indices);

figure(4); clf; hold on;
for jj = 1:3
    stim_levels = [0 1 2 3];
    for ii = 1:num_samples
        plot(stim_levels,squeeze(paramSweepResults.firing_rates_sorted(ii,jj,:)),'-.',...
            'color',[color_scheme(jj,:), 0.3],'linewidth',0.5);
    end
    
    mean_est = mean(squeeze(paramSweepResults.firing_rates_sorted(1:num_samples,jj,:)));
    h(jj) = plot(stim_levels,mean_est,'color',color_scheme(jj,:),'linewidth',4);
   
    set(gca,'fontsize',16)
    xticks([0 1 2 3])
    xticklabels({'None','Low','Med','High'})
    xlabel('Stimulus Strength')
    ylabel('Firing rate (Hz)')
end

legend([h(1) h(2) h(3)],{'PN','PV','SOM'})
box off

%% Plot 3f (distribution of successful recovery currents)
titles = {'PN','PV','SOM'};

figure(5); clf; 
% Set the edges for the histogram
histEdges = linspace(-5,5,10);
histEdges = [histEdges - diff(histEdges(1:2))/2, histEdges(end)+diff(histEdges(1:2))/2];
for ii = 1:3
    subplot(1,3,ii)
    histogram(paramSweepResults.recov_amounts_sorted(paramSweepResults.improved_gain_indices,ii), histEdges,...,
        'Normalization','pdf','Facecolor',color_scheme(ii,:),'Edgecolor',color_scheme(ii,:))
    
    set(gca,'fontsize',16)
    xlabel('I_{recov}')
    xticks([-4:2:4])
    title(titles{ii});
    box off
    
    numVals = length(paramSweepResults.improved_gain_indices);
    fprintf(strcat(titles{ii},' Mean: %.2f, STE: %.2f\n'), nanmean(paramSweepResults.recov_amounts_sorted(paramSweepResults.improved_gain_indices,ii)),...
        nanstd(paramSweepResults.recov_amounts_sorted(paramSweepResults.improved_gain_indices,ii))/sqrt(numVals));
end
