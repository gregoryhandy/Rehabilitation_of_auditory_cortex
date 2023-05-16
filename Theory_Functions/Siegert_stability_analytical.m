%%
% This is the analytical derivative of the Siegert function, as found in
% Helias, M., Tetzlaff T., and Diesmann, M. (2013; Appendix)
%
% Contains both the drates/dmu and drates/sigma terms
%%
function [W_stability] = Siegert_stability_analytical(firing_rates, mu_ext, var_ext,...
    tau_m,tau_r,V_th,V_r,J,Npop,J_sigma,alpha,tau_s,W)

m0 = J*firing_rates;
mu = (m0*tau_m + mu_ext');

var_0 = J_sigma*firing_rates;
sigma = sqrt(var_0*tau_m + var_ext');

alpha_Cee = zeros(Npop,1);
beta_Cee = zeros(Npop,1);
%% Find Cee
for j = 1:Npop
    % take the negative to convert 1+erf(x) to erfc below
    y_r_neg = -((V_r-mu(j))/sigma(j)+alpha/2*sqrt(tau_s/tau_m));
    y_th_neg = -((V_th-mu(j))/sigma(j)+alpha/2*sqrt(tau_s/tau_m));
    
    nu = firing_rates(j);
    
    front_factor = (nu*tau_m)^2*sqrt(pi)/sigma(j);
    term1 = exp(y_th_neg^2)*erfc(y_th_neg);
    term2 = exp(y_r_neg^2)*erfc(y_r_neg);
    alpha_Cee(j) = front_factor*(term1-term2);
    
    front_factor = sqrt(pi)*(tau_m*nu)^2/(2*sigma(j)^2);
    term1 = exp(y_th_neg^2)*erfc(y_th_neg)*(V_th-mu(j))/sigma(j);
    term2 = exp(y_r_neg^2)*erfc(y_r_neg)*(V_r-mu(j))/sigma(j);
    
    beta_Cee(j) = front_factor*(term1-term2);
end

% The first term comes from drates/dmu while the second is from
% drates/dsigma
W_stability = diag(alpha_Cee)*J+diag(beta_Cee)*(W.*J);

% P = inv(eye(params.Npop,params.Npop)-W_corr);



end

