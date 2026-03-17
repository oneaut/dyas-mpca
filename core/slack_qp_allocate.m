function [a_k, s_k, p_star] = slack_qp_allocate(tau_des, Phi_k, D1, P)
% SLACK_QP_ALLOCATE  Inner-loop reduced slack-QP allocation.
%
%   Solves the DyAS-reduced allocation problem (Eqs. 24-26 in paper):
%
%     min_{p, s}  p'*Hp*p + fp'*p + s'*W_tau*s
%     s.t.        B1*p + s = tau_des - B_phi*a_fix     (torque equality)
%                 D1*p in [a_lb - a_fix, a_ub - a_fix] (actuator bounds)
%
%   Inputs:
%     tau_des  - n x 1  desired joint torque from outer loop
%     Phi_k    - nw x nw diagonal deterioration matrix
%     D1       - nw x l  active subspace basis
%     P        - system parameter struct
%
%   Outputs:
%     a_k      - nw x 1  full actuator command (clipped to bounds)
%     s_k      - n x 1   slack (torque residual)
%     p_star   - l x 1   optimal reduced activation coefficients

n  = P.n;
nw = P.nw;
l  = size(D1, 2);

%% Compute complementary basis D2 via QR decomposition
[D_full, ~] = qr(D1);          % nw x nw orthonormal
D2 = D_full(:, l+1:end);       % nw x (nw-l)  inactive basis

%% Fixed inactive component
a_fix = D2 * (D2' * P.a_nom);  % nw x 1

%% Effective torque map and reduced quantities
B_phi = P.b * Phi_k;                  % n x nw
B1    = B_phi * D1;                   % n x l   reduced torque map
tau_rhs = tau_des - B_phi * a_fix;    % n x 1   adjusted torque target

%% QP cost matrices in reduced variable y = [p; s]
%   p block: Hp = D1'*R_u*D1 + eps*I
%   s block: W_tau (soft torque penalty)
Hp  = D1' * P.R_u * D1 + P.eps_Hp * eye(l);   % l x l
fp  = -2 * D1' * P.R_u * (P.a_nom - a_fix);    % l x 1

%   Full QP Hessian: block-diagonal [Hp, 0; 0, W_tau]
H_qp = blkdiag(Hp, P.W_tau);
f_qp = [fp; zeros(n, 1)];

%% Equality constraint:  B1*p + s = tau_rhs
%   [B1, I_n] * [p; s] = tau_rhs
Aeq = [B1, eye(n)];
beq = tau_rhs;

%% Inequality constraint: D1*p >= a_lb - a_fix  AND  D1*p <= a_ub - a_fix
%   Rewrite as: [-D1; D1] * p <= [-(a_lb-a_fix); (a_ub-a_fix)]
%   Pad with zeros for the slack variables
lb_a = P.a_lb - a_fix;   % nw x 1
ub_a = P.a_ub - a_fix;   % nw x 1

A_ineq = [-D1, zeros(nw, n);
           D1, zeros(nw, n)];
b_ineq = [-lb_a; ub_a];

%% Solve QP using quadprog (interior-point-convex)
opts = optimoptions('quadprog', ...
    'Algorithm',  'interior-point-convex', ...
    'Display',    'off', ...
    'MaxIterations', 200);

y0 = zeros(l + n, 1);   % initial guess

[y_star, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, ...
                                  Aeq, beq, [], [], y0, opts);

if exitflag <= 0
    % Fallback: use initial guess (slack absorbs all torque error)
    warning('slack_qp_allocate: quadprog exitflag = %d, using fallback.', exitflag);
    y_star = y0;
    y_star(l+1:end) = tau_rhs;   % put everything in slack
end

%% Extract solution
p_star = y_star(1:l);
s_k    = y_star(l+1:end);

%% Reconstruct full actuator vector and clip to bounds
a_k = D1 * p_star + a_fix;
a_k = max(P.a_lb, min(P.a_ub, a_k));

end
