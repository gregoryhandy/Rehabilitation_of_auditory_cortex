%%
% Make a raster plot 
% 1000 neurons are plotted, split between the population
%%
function plot_raster(times, tinds, Ntot, Npop, Ncells, pinds, plot_num )

% color_scheme = [0.8500 0.3250 0.0980;
%     0 0.4470 0.7410;
%     0.4660 0.6740 0.1880;
%     0.3 0.3 0.3];

% color_scheme = [81,128,68;
%     255,0,0;
%     221,124,1;
%     2,190,190]/255;

% color_scheme =[59, 57, 60; 164, 71, 105; 181, 117, 51]/255;

color_scheme =[59, 57, 60; 164, 71, 105;181, 117, 51; 107, 159, 165]/255;


%%

s = [times; tinds];
s=s(:,s(1,:)>0); % eliminate wasted space in s

% Time window over which to plot in ms
Tmin=0; Tmax=12000;

% Plot spikes of 1000 neurons
if Ntot > 100000
    plot_per_total = 1000/Ntot; % percent of neurons to plot
else
    plot_per_total = 1;
end
index_start = zeros(Npop,1);
index_end = zeros(Npop,1);
n_per_pop = zeros(Npop,1);
for i = 1:Npop
    % Plot cells from E population
%     n_per_pop(i) = Ncells(i)*plot_per_total;
    n_per_pop(i) = 100;
    
    index_start(i) = pinds(i);
    index_end(i) = pinds(i)+n_per_pop(i);
    
    Iplot=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:) > index_start(i) & s(2,:) <=index_end(i));
    
    % get the neuron indices
    neuroninds=s(2,Iplot);
    % map the neuron indices so that they lay next to each other
    if i > 1
        neuroninds = neuroninds - index_start(i) + sum(n_per_pop(1:(i-1)));
    end
    %neuroninds(neuroninds>Ne/2)=neuroninds(neuroninds>Ne/2)-Ne/2+nplot/2;
    
    if ~isempty(neuroninds)
        s_raster = double(s(1,Iplot))/1000;
        h(i) = plot(s_raster,neuroninds,'.','MarkerSize',6,'color', color_scheme(i,:));
        hold on;
    else
        % plot a point out of the axis so that the legend still works
        h(i) = plot(-10,-10,'.','MarkerSize',6,'color', color_scheme(i,:));
    end
end

xlabel('Time (sec)')
yticklabels([])
set(gca,'fontsize',16)
    
if plot_num == 4
    if Npop == 4
        [lgd,icons] = legend([h(4), h(3), h(2), h(1)],{'VIP','SOM','PV','PN'},'fontsize',16);
    elseif Npop == 3
        if isempty(neuroninds)
            h(i) = plot(-1,-1,'.','MarkerSize',6,'color', color_scheme(i,:));
        end
        [lgd,icons] = legend([h(3), h(2), h(1)],{'SOM','PV','PN'},'fontsize',16);
    end

    if Npop>=3
        % Create a legend with 3 entries
        % Find the 'line' objects
        icons = findobj(icons,'Type','line');
        % % Find lines that use a marker
        icons = findobj(icons,'Marker','none','-xor');
        % % Resize the marker in the legend
        set(icons,'MarkerSize',20);
        % lgd.FontSize = 100;
    end
end

ylim([0 sum(n_per_pop)])

end



