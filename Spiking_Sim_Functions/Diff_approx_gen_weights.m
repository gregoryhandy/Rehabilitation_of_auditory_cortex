%%
% Create the connection matrix used in the diffusion approximation
% Note: There are no feedforward connections, since these are approximated
% in the main function
%%
function [wind,wipost,wstr] = Diff_approx_gen_weights(Ncells,p0,J)

%% set up recurrent weight matrix
Ntot = sum(Ncells);
Npop = length(Ncells);

Maxw = round(Ntot*Ntot*0.3); % maximum number of weights in the weight matrix
wind = zeros(Ntot+1,1); % column of w corresponding to the start of the ith neuron's projections                
	
wipost = zeros(Maxw,1);
wstr = zeros(Maxw,1);

syncount = 1;
% loop through the populations
for pp = 1:Npop
    
    % loop through each neuron
    for cc=1:Ncells(pp)
        
        starting_Index = sum(Ncells(1:(pp-1)));
        wind(cc + starting_Index) = syncount;
        
        % find which neurons are connected to neuron cc
        for qq = 1:Npop
        
            % probability of a connection
            prob = p0(pp,qq);
            
            % Method 1: flip weighted coins!            
            % iconns = find(rand(Ncells(qq),1) < prob) + sum(Ncells(1:(qq-1)));
            
            % Method 2: fix the total number of connecctions to be Ncells(qq)*prob
            iconns = randperm(Ncells(qq),round(Ncells(qq)*prob))+sum(Ncells(1:(qq-1)));
             
            % record the connections
            wipost(syncount:(syncount+length(iconns)-1)) = iconns;
            wstr(syncount:(syncount+length(iconns)-1)) = J(pp,qq);
            
            syncount = syncount+length(iconns);
        end
    end
end

wind(Ntot+1)= syncount-1;
% reshape these vectors to get rid of unnecessary zeros
wipost = wipost(1:(syncount-1));
wstr = wstr(1:(syncount-1));

end