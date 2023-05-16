%%
% Recreates SFig 5 for Kumar et al., 2023
%
% Compares the firing rates for the pre- and post-damage cases for the
% default circuit configuration and the network with no feedforward inputs
% onto VIPs
% The feedforward inputs are essentially "replace" by manually modeling
% short-term facilitation of excitatory inputs onto VIP neuron
%%
clear; close all; clc;

% Load the datasets
load('./Sim_Data/SFig5Data')

%% SFig 5a (default circuit configuration)
color_scheme =[59, 57, 60; 164, 71, 105;181, 117, 51; 107, 159, 165]/255;
f=figure(1); clf; f.Position(3:4) = [1200, 400];
subplot(1,3,1); hold on;
for jj = 1:postDefault.params.Npop
    hLeg(jj) = plot(1:1:4,preDefault.firing_rates_sol(jj,:),'.-','markersize',30,...
        'color',color_scheme(jj,:),'linewidth',2);  
    plot(1:1:4,postDefault.firing_rates_sol(jj,:),'s--','markersize',10,...
        'MarkerFaceColor',color_scheme(jj,:),'linewidth',2,'color',color_scheme(jj,:));  
end
set(gca,'fontsize',16)
legend(hLeg, {'PN','PV','SOM','VIP'})
xlabel('Stimulus Strength')
ylabel('Firing Rate (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})
box off
ylim([0 40])
title('Default Circuit')

%% SFig 5b (replace the ffwd input on VIP with short-term exc. facilitation)
subplot(1,3,2); hold on;
for jj = 1:postDefault.params.Npop
    gLeg(jj) = plot([1:1:4],preAdj.firing_rates_sol(jj,:),'.-','markersize',30,...
        'color',color_scheme(jj,:),'linewidth',2);  
    hLeg(jj) = plot([1:1:4],postAdj.firing_rates_sol(jj,:),'s--','markersize',10,...
        'MarkerFaceColor',color_scheme(jj,:),'linewidth',2,'color',color_scheme(jj,:));
end
set(gca,'fontsize',16)
legend([gLeg(1) hLeg(1)], {'Pre','Post'})
xlabel('Stimulus Strength')
ylabel('Firing Rate (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})
ylim([0 40])
box off
title('No Ffwd Input onto VIPs')

%% SFig 5c (compare panels a and b) 
subplot(1,3,3); hold on;
for jj = 1:postDefault.params.Npop
    plot(1:1:4,abs(postDefault.firing_rates_sol(jj,:)-postAdj.firing_rates_sol(jj,:)),...
        's--','markersize',10,'MarkerFaceColor',color_scheme(jj,:),'linewidth',2,'color',color_scheme(jj,:));
    plot(1:1:4,abs(preDefault.firing_rates_sol(jj,:)-preAdj.firing_rates_sol(jj,:)),'.-',...
        'markersize',30,'color',color_scheme(jj,:),'linewidth',2);
end
ylim([0 40])
set(gca,'fontsize',16)
xlabel('Stimulus Strength')
ylabel('Abs Firing Rate Difference (Hz)')
xticks([1 2 3 4])
xticklabels({'None','Low','Med','High'})

title(sprintf('Difference between the \n two configurations'))