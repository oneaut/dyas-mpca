function [q_ref, dq_ref, ddq_ref] = generate_sts_reference(P, n_steps, dt)
% GENERATE_STS_REFERENCE  Minimum-jerk reference for 20-cycle STS session.
%
%   Generates a sequence of minimum-jerk trajectories for the STS cycle:
%     Phase 1 (0-3s):  sit-to-stand  q0 -> q1
%     Phase 2 (3-4s):  standing hold q1 -> q1
%     Phase 3 (4-7s):  stand-to-sit  q1 -> q0
%   Repeated for P.N_cycles cycles.
%
%   Reference: Flash & Hogan, J. Neurosci., 1985.

n_cycle = round(P.T_cycle / dt);
n_rise  = round(P.t_rise  / dt);
n_hold  = round(P.t_hold  / dt);
n_sit   = round(P.t_sit   / dt);

% Pre-allocate
q_ref   = zeros(P.n, n_steps);
dq_ref  = zeros(P.n, n_steps);
ddq_ref = zeros(P.n, n_steps);

% Build one cycle
for phase = 1:3
    switch phase
        case 1   % sit-to-stand
            q_start = P.q0;  q_end = P.q1;  T = P.t_rise;  ns = n_rise;
        case 2   % hold
            q_start = P.q1;  q_end = P.q1;  T = P.t_hold;  ns = n_hold;
        case 3   % stand-to-sit
            q_start = P.q1;  q_end = P.q0;  T = P.t_sit;   ns = n_sit;
    end

    if phase == 1
        offset = 0;
    elseif phase == 2
        offset = n_rise;
    else
        offset = n_rise + n_hold;
    end

    tau_vec = linspace(0, T, ns);
    for k = 1:ns
        tau_n = tau_vec(k) / T;  % normalised time in [0,1]
        % Minimum-jerk position, velocity, acceleration
        s   = 10*tau_n^3 - 15*tau_n^4 + 6*tau_n^5;
        ds  = (30*tau_n^2 - 60*tau_n^3 + 30*tau_n^4) / T;
        dds = (60*tau_n   - 180*tau_n^2 + 120*tau_n^3) / T^2;

        q_ref(:, offset+k)   = q_start + s   * (q_end - q_start);
        dq_ref(:, offset+k)  = ds  * (q_end - q_start);
        ddq_ref(:, offset+k) = dds * (q_end - q_start);
    end
end

% One-cycle template
q_one   = q_ref(:,   1:n_cycle);
dq_one  = dq_ref(:,  1:n_cycle);
ddq_one = ddq_ref(:, 1:n_cycle);

% Tile for all cycles
for c = 1:P.N_cycles
    idx = (c-1)*n_cycle + (1:n_cycle);
    if idx(end) > n_steps
        idx = idx(1):n_steps;
    end
    q_ref(:,   idx) = q_one(:,   1:numel(idx));
    dq_ref(:,  idx) = dq_one(:,  1:numel(idx));
    ddq_ref(:, idx) = ddq_one(:, 1:numel(idx));
end

end
