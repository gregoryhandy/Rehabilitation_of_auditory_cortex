function [rates_ave] = spikeAnalysisFn(tinds,times,params)

tBurn = 1000;
rates_ave = zeros(params.Npop,1);
for jj = 1:params.Npop
    % indicator for the current population
    popData = ismember(tinds,params.pinds(jj):params.pinds(jj+1)-1);
    % indicator for the right time window
    timeData = times>=tBurn;
    
    % Add up total spikes and average by num cells, time and convert to Hz
    rates_ave(jj) = sum(popData & timeData)/(params.Ncells(jj)*(params.T-tBurn))*1e3;
end

end

