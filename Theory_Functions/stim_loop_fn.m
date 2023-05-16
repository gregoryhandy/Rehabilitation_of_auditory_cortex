%%
% Loops over the strengths of the stimuli considered
% Default is fours values (none, low, med, and high)
%
% Returns the firing rates (of the mean field model) and the eigenvalues to
% assess the stability of the system
%%
function [firing_rates_sol,nan_warning,warning_notice, eig_values, d_min, real_eig_max]...
    = stim_loop_fn(params,bg_damage,stim_damage,recov_amounts)

%% Preallocate
firing_rates_sol = zeros(params.Npop,params.num_stims);
warning_notice = zeros(params.num_stims,1);
nan_warning = zeros(params.num_stims,1);

d_min = zeros(params.num_stims,1);
real_eig_max = zeros(params.num_stims,1);
eig_values = zeros(params.Npop,params.num_stims);
omega = [0:0.1:1000];

%% Loop over all stimulus
for ii  = 1:params.num_stims
    
    %% Estimate the steady state firing rate
    mu_ext = params.mu_stim(ii,:).*stim_damage + ...
        params.mu_bg.*bg_damage+recov_amounts;
    var_ext = params.sigma_stim(ii,:).^2.*stim_damage+...
        params.sigma_bg.^2.*bg_damage+params.sigma_fixed.^2;
    
    [firing_rates, nan_warning(ii), warning_notice(ii)] = ...
        theory_firing_rate_fn(...
        mu_ext, var_ext,params.tau_m_theory,params.tau_ref_theory,...
        params.V_th,params.V_r, params.J_theory,...
        params.Npop,params.J_sigma_theory,params.alpha, ...
        params.tau_s_theory);
      
    % Correct for numerical error (for the stability measurement)
    firing_rates(firing_rates<0.0001) = 0;
    firing_rates_sol(:,ii) = firing_rates;
    
    %% Find the stability of the fixed point
    [W_stability] = Siegert_stability_analytical(firing_rates_sol(:,ii), ...
        mu_ext, var_ext,...
        params.tau_m_theory,params.tau_ref_theory,...
        params.V_th,params.V_r, params.J_theory,...
        params.Npop,params.J_sigma_theory,params.alpha, ...
        params.tau_s_theory,params.W);
    
    eig_values(:,ii) = eig(W_stability-eye(params.Npop,params.Npop));
    lambda_fn = eig_values(:,ii)./(1+1i*omega);
    
    d_min(ii) = min(min(abs(1-lambda_fn)));
    real_eig_max(ii) = max(real(eig_values(:,ii)));
end

end

