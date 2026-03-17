function P = sts_params()
% STS_PARAMS  All parameters for the 3-DOF STS hybrid neuroprosthesis.
%
%   P = sts_params() returns a struct P containing all system, controller,
%   and DyAS parameters used in:
%
%   M. Singh and N. Sharma, "Stability-Guaranteed Dynamic Active Subspace
%   Torque Allocation for Over-Actuated Hybrid Neuroprostheses,"
%   IEEE Control Systems Letters, 2025.

% -------------------------------------------------------------------------
%  Dimensions
% -------------------------------------------------------------------------
P.n  = 3;   % DOF (ankle, knee, hip)
P.nw = 8;   % Actuators (3 motors + 5 FES channels)
P.nm = 3;   % Number of motors
P.nf = 5;   % Number of FES channels

% -------------------------------------------------------------------------
%  Anthropometrics  (Ashby et al., J. Rehabil. Res. Dev., 2004)
% -------------------------------------------------------------------------
P.m_seg = [1.2, 3.7, 8.6, 40.0];      % Segment masses [kg]: foot/shank/thigh/HAT
P.l_seg = [0.26, 0.43, 0.42, 0.55];   % Segment lengths [m]

% -------------------------------------------------------------------------
%  FES channel parameters  (order: TA, SOL, QUAD, HAMS, GLUT)
% -------------------------------------------------------------------------
P.F_max  = [800, 1200, 1500, 1000, 1400];  % Max isometric forces [N]

% Moment arms [m]
P.r_FES = [ 0.04, -0.05,  0.00,  0.00,  0.00;   % Ankle
            0.00,  0.00,  0.04, -0.04,  0.00;   % Knee
            0.00,  0.00,  0.00, -0.05, -0.06];  % Hip

% Torque effectiveness matrix  b = [I_{3x3} | R_FES * diag(F_max)]
R_FES_scaled = P.r_FES * diag(P.F_max);
P.b = [eye(P.nm), R_FES_scaled];   % 3 x 8

% Baseline fatigue and recovery rates [s^-1]  (TA, SOL, QUAD, HAMS, GLUT)
P.wf0 = [0.018, 0.015, 0.025, 0.020, 0.022];
P.wr0 = [0.009, 0.008, 0.013, 0.010, 0.011];

P.phi_min = 0.20 * ones(1, P.nf);   % Minimum active fraction
P.phi_max = 1.00 * ones(1, P.nf);   % Maximum active fraction

% -------------------------------------------------------------------------
%  Motor and actuator bounds
% -------------------------------------------------------------------------
P.tau_max_motor = [200, 250, 200];  % Peak motor torques [Nm]: ankle/knee/hip

P.a_lb  = [-1, -1, -1,  0,  0,  0,  0,  0]';   % 8 x 1
P.a_ub  = [ 1,  1,  1,  1,  1,  1,  1,  1]';   % 8 x 1
P.a_nom = zeros(P.nw, 1);                        % Nominal (rest) input

% -------------------------------------------------------------------------
%  Passive joint damping  (diagonal, used in CT-PD and dynamics)
% -------------------------------------------------------------------------
P.D_damp = diag([5, 8, 6]);   % [Nm/(rad/s)]: ankle, knee, hip

% -------------------------------------------------------------------------
%  STS motion parameters
% -------------------------------------------------------------------------
P.q0 = [5,  40, -35]' * pi/180;   % Initial seated posture [rad]
P.q1 = [3,   5,   0]' * pi/180;   % Standing target posture [rad]

P.t_rise   = 3.0;   % Sit-to-stand duration [s]
P.t_hold   = 1.0;   % Standing hold duration [s]
P.t_sit    = 3.0;   % Stand-to-sit duration [s]
P.T_cycle  = P.t_rise + P.t_hold + P.t_sit;   % 7 s per cycle
P.N_cycles = 20;                                % Cycles per session
P.T_session = P.N_cycles * P.T_cycle;           % 140 s total

% -------------------------------------------------------------------------
%  Controller gains
% -------------------------------------------------------------------------
P.Kp = diag([500, 450, 400]);    % Proportional gain
P.Kd = diag([ 45,  42,  38]);    % Derivative gain

P.W_tau  = diag([10, 10, 10]);                           % Torque tracking weight  3x3
P.R_u    = diag([0.05, 0.05, 0.05, 2, 2, 2, 2, 2]);    % Effort weight           8x8
P.eps_Hp = 1e-6;                                          % Hessian regularisation

% -------------------------------------------------------------------------
%  Simulation settings
% -------------------------------------------------------------------------
P.dt           = 0.01;   % Control timestep [s] = 100 Hz
P.rk4_substeps = 4;      % RK4 substeps per control step

end
