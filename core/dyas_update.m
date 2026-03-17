function [D1, l_star, C_new, Pi_new] = dyas_update(Phi_k, tau_des, C_prev, Pi_prev, P, D)
% DYAS_UPDATE  Online DyAS identification using analytical gradient.
%   Algorithm 1 from Singh & Sharma, IEEE CSL 2025.
%
%   Inputs:
%     Phi_k    - nw x nw diagonal deterioration matrix at step k
%     tau_des  - n x 1 desired joint torque
%     C_prev   - nw x nw previous GCM
%     Pi_prev  - nw x nw previous eigenvector propagation matrix
%     P        - system parameter struct (from sts_params)
%     D        - DyAS parameter struct  (from dyas_params)
%
%   Outputs:
%     D1       - nw x l_star active subspace basis
%     l_star   - selected active dimension
%     C_new    - nw x nw updated GCM
%     Pi_new   - nw x nw updated propagation matrix

nw = P.nw;
Ms = D.M_s;

%% Step 1: Effective torque map
B_phi = P.b * Phi_k;                           % n x nw

%% Step 2: Cost Hessian and gradient offset (once per update)
%   H = B_phi' * W_tau * B_phi + R_u          (nw x nw)
%   r = B_phi' * W_tau * tau_des + R_u * a_nom (nw x 1)
H = B_phi' * P.W_tau * B_phi + P.R_u;          % nw x nw, positive definite
r = B_phi' * P.W_tau * tau_des + P.R_u * P.a_nom;  % nw x 1

%% Step 3: Sample random actuator points
U_s = P.a_lb + (P.a_ub - P.a_lb) .* rand(nw, Ms);  % nw x Ms

%% Step 4: Analytical gradient matrix (NO plant evaluations)
%   G = 2*H*U_s - 2*r*1'       (nw x Ms)
G = 2 * H * U_s - 2 * r * ones(1, Ms);

%% Step 5: Batch GCM estimate
C_hat = (1 / Ms) * G * G';        % nw x nw

%% Step 6: Update eigenvector propagation matrix Pi
%   Least-squares estimate:  Pi * D_prev = D_new
%   Updated incrementally to track eigenvector rotation.
[D_prev, ~] = eig(C_prev, 'vector');
D_prev = fliplr(D_prev);   % descending eigenvalue order

[D_hat, ~]  = eig(C_hat, 'vector');
D_hat = fliplr(D_hat);

% Incremental update of Pi
Pi_new = Pi_prev;
for j = 1:nw
    d_prev_j = Pi_prev * D_prev(:,j);
    nrm = norm(d_prev_j);
    if nrm > 1e-10
        Pi_new = Pi_new + (D_hat(:,j) - d_prev_j) * D_prev(:,j)' / (nrm^2);
    end
end

%% Step 7: EMA update of GCM
C_new = D.mu * C_hat + (1 - D.mu) * Pi_new * C_prev * Pi_new';

% Symmetrise to prevent numerical drift
C_new = 0.5 * (C_new + C_new');

%% Step 8: Eigendecomposition of updated GCM
[V, lam_vec] = eig(C_new, 'vector');

% Sort descending
[lam_sorted, idx] = sort(lam_vec, 'descend');
V_sorted = V(:, idx);

% Clip small negative eigenvalues from numerical errors
lam_sorted = max(lam_sorted, 0);

%% Step 9: Select active dimension via variance threshold
cumvar  = cumsum(lam_sorted) / sum(lam_sorted);
l_star  = find(cumvar >= D.eta, 1, 'first');

% Safety clamp
l_star = max(D.l_min, min(D.l_max, l_star));

%% Step 10: Return active basis
D1 = V_sorted(:, 1:l_star);     % nw x l_star

end
