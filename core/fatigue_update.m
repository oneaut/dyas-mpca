function phi_new = fatigue_update(phi, a, wf, wr, P, dt)
% FATIGUE_UPDATE  Euler integration of FES fatigue dynamics (Eq. 4 in paper).
%
%   dphi_i/dt = wf_i*(phi_min_i - phi_i)*a_tilde_i
%             + wr_i*(1 - phi_i)*(1 - a_tilde_i)
%
%   Inputs:
%     phi  - nw x 1  current active fractions  (motors: fixed at 1)
%     a    - nw x 1  actuator activations
%     wf   - nw x 1  fatigue rate vector  (0 for motors)
%     wr   - nw x 1  recovery rate vector (0 for motors)
%     P    - system parameter struct
%     dt   - timestep [s]
%
%   Output:
%     phi_new - nw x 1  updated active fractions

nw = P.nw;

% Normalised activation
a_tilde = a ./ P.a_ub;                 % nw x 1
a_tilde = max(0, min(1, a_tilde));      % clip to [0,1]

% Build phi_min vector (0 for motors, phi_min for FES channels)
phi_min_vec = [zeros(P.nm, 1); P.phi_min(:)];

% Fatigue dynamics
dphi = wf .* (phi_min_vec - phi) .* a_tilde ...
     + wr .* (1 - phi)           .* (1 - a_tilde);

% Euler integration
phi_new = phi + dt * dphi;

% Clamp to valid range
phi_min_full = [ones(P.nm, 1); P.phi_min(:)];   % motors always 1
phi_max_full = ones(nw, 1);

phi_new(1:P.nm) = 1.0;   % motors do not fatigue
phi_new = max(phi_min_full, min(phi_max_full, phi_new));

end
