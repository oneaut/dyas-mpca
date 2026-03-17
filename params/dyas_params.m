function D = dyas_params()
% DYAS_PARAMS  DyAS hyperparameters used in Singh & Sharma, IEEE CSL 2025.
%
%   D = dyas_params() returns a struct D with all DyAS-specific settings.

D.M_s    = 20;     % Number of gradient samples per GCM update
D.mu     = 0.30;   % EMA smoothing weight for GCM update
D.eta    = 0.95;   % Variance threshold for active dimension selection
D.dt_dyas = 0.10;  % DyAS update period [s]  (every 10 control steps at 100 Hz)

% Active subspace dimension bounds (safety clamp)
D.l_min  = 1;
D.l_max  = 8;   % = n_w; fallback to full if GCM is ill-conditioned

end
