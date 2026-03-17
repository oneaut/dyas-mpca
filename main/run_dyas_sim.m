function run_dyas_sim(varargin)
% RUN_DYAS_SIM  Run the full DyAS-MPCA simulation for one or all fatigue scenarios.
%
%   run_dyas_sim()
%       Runs all three scenarios (A, B, C) and saves results to results/
%
%   run_dyas_sim('scenario', 'B_rapid')
%       Runs a single named scenario.
%
%   Scenarios: 'A_mild'  |  'B_rapid'  |  'C_prefatigued'
%
%   Reference:
%   M. Singh and N. Sharma, "Stability-Guaranteed Dynamic Active Subspace
%   Torque Allocation for Over-Actuated Hybrid Neuroprostheses,"
%   IEEE Control Systems Letters, 2025.

%% Parse inputs
p = inputParser;
addParameter(p, 'scenario', 'all');
parse(p, varargin{:});
scenario_arg = p.Results.scenario;

%% Add all subfolders to path
rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(rootdir));

%% Choose which scenarios to run
if strcmp(scenario_arg, 'all')
    scenario_list = {'A_mild', 'B_rapid', 'C_prefatigued'};
else
    scenario_list = {scenario_arg};
end

%% Load fixed system parameters
P = sts_params();
D = dyas_params();

%% Run each scenario
for s = 1:numel(scenario_list)
    sc_name = scenario_list{s};
    S = fatigue_scenarios(sc_name);

    fprintf('\n=== Running scenario: %s ===\n', S.label);

    % Run DyAS allocator
    fprintf('  DyAS allocator ...\n');
    results_dyas = run_one_session(P, D, S, 'DyAS');

    % Run Full allocator
    fprintf('  Full allocator  ...\n');
    results_full = run_one_session(P, D, S, 'Full');

    % Save
    outfile = fullfile(rootdir, 'results', ['results_', sc_name, '.mat']);
    save(outfile, 'results_dyas', 'results_full', 'P', 'D', 'S');
    fprintf('  Saved: %s\n', outfile);
end

fprintf('\nAll scenarios complete.\n');
fprintf('Run compare_results() to generate figures.\n');
end


% =========================================================================
function R = run_one_session(P, D, S, mode)
% Run one complete 20-cycle STS session.

dt      = P.dt;
n_steps = round(P.T_session / dt);

%% Initialise state
q    = P.q0;
dq   = zeros(P.n, 1);
phi  = [ones(P.nm, 1); S.phi0(:)];   % [motors; FES channels]

%% Fatigue rates for this scenario
wf = [zeros(P.nm, 1); (S.wf_mult * P.wf0)'];
wr = [zeros(P.nm, 1); (S.wr_mult * P.wr0)'];

%% Initialise DyAS state
C_gcm = eye(P.nw);          % GCM
Pi    = eye(P.nw);           % Eigenvector propagation matrix
D1    = eye(P.nw);           % Active subspace basis (full at start)
l_star = P.nw;
t_last_dyas = -inf;

%% Pre-generate minimum-jerk reference for all cycles
[q_ref, dq_ref, ddq_ref] = generate_sts_reference(P, n_steps, dt);

%% Storage
store_q    = zeros(P.n,  n_steps);
store_dq   = zeros(P.n,  n_steps);
store_phi  = zeros(P.nw, n_steps);
store_a    = zeros(P.nw, n_steps);
store_tau  = zeros(P.n,  n_steps);
store_l    = zeros(1,    n_steps);
store_s    = zeros(P.n,  n_steps);
store_rmse = zeros(P.n,  n_steps);

%% Main simulation loop
for k = 1:n_steps
    t = (k-1) * dt;
    Phi_k = diag(phi);

    %% 1. Outer loop: computed-torque PD
    tau_des = ctpd_outer_loop(q, dq, q_ref(:,k), dq_ref(:,k), ...
                              ddq_ref(:,k), P, Phi_k);

    %% 2. DyAS update (if due)
    if strcmp(mode, 'DyAS') && (t - t_last_dyas) >= D.dt_dyas
        [D1, l_star, C_gcm, Pi] = dyas_update(Phi_k, tau_des, ...
                                               C_gcm, Pi, P, D);
        t_last_dyas = t;
    end

    %% 3. Inner loop: allocate torque
    if strcmp(mode, 'DyAS')
        [a_k, s_k] = slack_qp_allocate(tau_des, Phi_k, D1, P);
    else
        % Full allocation: D1 = I, l* = nw
        D1_full = eye(P.nw);
        [a_k, s_k] = slack_qp_allocate(tau_des, Phi_k, D1_full, P);
        l_star = P.nw;
    end

    %% 4. Apply torque and integrate plant
    tau_applied = P.b * Phi_k * a_k;
    [q_new, dq_new] = sts_dynamics(q, dq, tau_applied, P, dt);

    %% 5. Update fatigue
    phi = fatigue_update(phi, a_k, wf, wr, P, dt);

    %% 6. Store
    store_q(:,k)    = q;
    store_dq(:,k)   = dq;
    store_phi(:,k)  = phi;
    store_a(:,k)    = a_k;
    store_tau(:,k)  = tau_applied;
    store_l(k)      = l_star;
    store_s(:,k)    = s_k;
    store_rmse(:,k) = (q - q_ref(:,k)).^2;

    %% Advance state
    q  = q_new;
    dq = dq_new;
end

%% Package results
R.q      = store_q;
R.dq     = store_dq;
R.phi    = store_phi;
R.a      = store_a;
R.tau    = store_tau;
R.l_star = store_l;
R.s      = store_s;
R.rmse   = sqrt(store_rmse) * (180/pi);   % convert to degrees
R.q_ref  = q_ref;
R.t      = (0:n_steps-1) * dt;
R.mode   = mode;

% Per-cycle summary
n_steps_cycle = round(P.T_cycle / dt);
R.phi_per_cycle  = zeros(P.nw, P.N_cycles);
R.rmse_per_cycle = zeros(P.n,  P.N_cycles);
for c = 1:P.N_cycles
    idx = (c-1)*n_steps_cycle + (1:n_steps_cycle);
    R.phi_per_cycle(:,c)  = mean(store_phi(:,idx), 2);
    R.rmse_per_cycle(:,c) = mean(sqrt(store_rmse(:,idx)), 2) * (180/pi);
end
end
