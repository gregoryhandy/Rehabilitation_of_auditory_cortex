%%
% Population steady state calculation using adjusted Eqn. 4.33 
% from de la Rocha, 2007 (converts the integral to use erfc)
%%
function [firing_rates,nan_warning,warning_notice] = ...
    theory_firing_rate_fn(m_ext, var_ext,tau_m,tau_r,V_th,V_r,J,Npop,...
    J_sigma,alpha,tau_s)


%% Find the steady state firing rates
max_change = 10;
tt = 1;

nan_warning = 0;
max_counter = 10000;

result = zeros(Npop,1);
dt = 0.01;

firing_rates_time = zeros(Npop,max_counter);

%% Initial firing rates
firing_rates_time(:,1) = 5*ones(Npop,1);

%% Loop over time to find the fixed point
while max_change > 1e-5 && tt < max_counter
    
    m0 = J*firing_rates_time(:,tt);
    mu = (m0*tau_m + m_ext');
   
    var_0 = J_sigma*firing_rates_time(:,tt);
    sigma = sqrt(var_0*tau_m + var_ext');
    
    y_th = -((V_th-mu)./sigma+alpha/2*sqrt(tau_s/tau_m));
    
    y_r = -((V_r-mu)./sigma+alpha/2*sqrt(tau_s/tau_m));
    
    nan_test_y_r = sum(isnan(y_r));
    nan_test_y_th = sum(isnan(y_th));
    if nan_test_y_r > 0 || nan_test_y_th > 0
       nan_warning = 1;
       break;
    end
    
    % loop through the populations
    for j = 1:Npop
        temp_x = [y_th(j):0.001:y_r(j)];
        result(j) = 1/(tau_r+trapz(temp_x,exp(temp_x.^2).*(erfc(temp_x)))*tau_m*sqrt(pi));
    end
    delta_y = -firing_rates_time(:,tt) + result;
    max_change = max(abs(delta_y));
%     max_change = max(abs(delta_y*dt));
    firing_rates_time(:,tt+1) = firing_rates_time(:,tt) + delta_y*dt;
    

    % start to worry about oscillations if the system hasn't converged by
    % this point
%     if tt*dt > 100
%         if firing_rates_time(1,tt+1) > max_fr
%             max_fr = firing_rates_time(1,tt+1);
%         end
%         
%         if firing_rates_time(1,tt+1) < min_fr
%             min_fr = firing_rates_time(1,tt+1);
%         end
%     end
    
    % if this is the case, most likely the point is unstable and we are
    % oscillating!
%     if (max_fr - min_fr) > 1
%         tt = max_counter;
%         break;
%     end
    
    % take all rates to be zero if this is the case
    % speeds up the algorithm
    if max(firing_rates_time(:,tt+1) ) < 1e-2
        firing_rates_time(:,tt)  = 0;
        break;
    end
    
    tt = tt+1;
end

if tt == max_counter
    warning_notice = max_change;
else
    warning_notice = 0;
end

firing_rates = firing_rates_time(:,tt);

% figure(99)
% hold off
% plot(firing_rates_time(1,:))
% hold on
% plot(firing_rates_time(2,:))
% plot(firing_rates_time(3,:))

end

