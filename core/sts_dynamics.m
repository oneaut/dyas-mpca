function [q_new, dq_new] = sts_dynamics(q, dq, tau, P, dt)
% STS_DYNAMICS  Integrate 3-DOF STS hybrid neuroprosthesis dynamics one step.
%
%   Applies RK4 with P.rk4_substeps substeps per control step.
%
%   Dynamics (Eq. 1 in paper):
%     M(q)*ddq + C(q,dq)*dq + G(q) + D*dq = tau
%
%   Inputs:
%     q    - n x 1  joint angles  [rad]
%     dq   - n x 1  joint velocities [rad/s]
%     tau  - n x 1  applied joint torques [Nm]
%     P    - system parameter struct
%     dt   - control timestep [s]
%
%   Outputs:
%     q_new   - n x 1  updated joint angles
%     dq_new  - n x 1  updated joint velocities

h  = dt / P.rk4_substeps;   % RK4 substep size

x = [q; dq];   % state vector: [angles; velocities]

for i = 1:P.rk4_substeps
    k1 = sts_xdot(x,           tau, P);
    k2 = sts_xdot(x + h/2*k1,  tau, P);
    k3 = sts_xdot(x + h/2*k2,  tau, P);
    k4 = sts_xdot(x + h*k3,    tau, P);
    x  = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

q_new  = x(1:P.n);
dq_new = x(P.n+1:end);

end


function xdot = sts_xdot(x, tau, P)
% State derivative: xdot = [dq; ddq]

q  = x(1:P.n);
dq = x(P.n+1:end);

M = sts_inertia(q, P);
C = sts_coriolis(q, dq, P);
G = sts_gravity(q, P);

% Equations of motion:  ddq = M^{-1}(tau - C*dq - G - D*dq)
ddq = M \ (tau - (C + P.D_damp) * dq - G);

xdot = [dq; ddq];
end


function M = sts_inertia(q, P)
% 3x3 inertia matrix for 4-segment sagittal STS model.
% Segments: foot(1), shank(2), thigh(3), HAT(4)
% Joints:   ankle(1), knee(2), hip(3)

m  = P.m_seg;
l  = P.l_seg;
lc = l / 2;    % centre-of-mass at segment midpoint (simplification)

% Segment moments of inertia about their CoM (thin rod approximation)
I = (1/12) * m .* l.^2;

% Absolute angles (each joint relative to vertical)
theta1 = q(1);   % ankle
theta2 = q(2);   % knee
theta3 = q(3);   % hip

% Using a simplified 3-DOF lumped parameter inertia
% (full expression omitted for brevity; replace with patient-specific model)
m2 = m(2); m3 = m(3); m4 = m(4);
l2 = l(2); l3 = l(3);
lc2 = lc(2); lc3 = lc(3);

I2 = I(2); I3 = I(3); I4 = I(4);

% Inertia contributions (simplified 3-DOF)
M11 = (m2+m3+m4)*l2^2 + (m3+m4)*l3^2 + I2 + I3 + I4 ...
      + 2*(m3+m4)*l2*l3*cos(theta2) ...
      + 2*m4*l2*lc3*cos(theta2) + 2*m4*l3*lc3*cos(theta3);

M12 = (m3+m4)*l3^2 + I3 + I4 ...
      + (m3+m4)*l2*l3*cos(theta2) + m4*l3*lc3*cos(theta3);

M13 = m4*lc3^2 + I4;

M22 = (m3+m4)*l3^2 + I3 + I4 + 2*m4*l3*lc3*cos(theta3);
M23 = m4*lc3^2 + I4 + m4*l3*lc3*cos(theta3);
M33 = m4*lc3^2 + I4;

M = [M11, M12, M13;
     M12, M22, M23;
     M13, M23, M33];

% Ensure positive definiteness
M = 0.5*(M + M') + 1e-6*eye(3);
end


function C = sts_coriolis(q, dq, P)
% 3x3 Coriolis/centripetal matrix (simplified)

l  = P.l_seg;
m  = P.m_seg;

theta2 = q(2);
theta3 = q(3);
dq1 = dq(1); dq2 = dq(2); dq3 = dq(3);

m3 = m(3); m4 = m(4);
l2 = l(2); l3 = l(3);

h2 = -(m3+m4)*l2*l3*sin(theta2);
h3 = -m4*l2*l3*sin(theta3);   % simplified

C = [h2*(2*dq2+dq3)+h3*dq3,  h2*dq1,           h3*dq1;
     -h2*dq1,                  0,                 h3*(dq1+dq2);
     -(h3*(2*dq1+dq2)),        -h3*dq2,           0];
end


function G = sts_gravity(q, P)
% 3x1 gravity vector

m  = P.m_seg;
l  = P.l_seg;
lc = l/2;
g  = 9.81;

m2 = m(2); m3 = m(3); m4 = m(4);
l2 = l(2); l3 = l(3);
lc3 = lc(3);

theta1 = q(1);
theta2 = q(2);
theta3 = q(3);

% Gravity contributions projected onto joint torques
G1 = (m2+m3+m4)*g*l2*cos(theta1) + (m3+m4)*g*l3*cos(theta1+theta2) ...
     + m4*g*lc3*cos(theta1+theta2+theta3);

G2 = (m3+m4)*g*l3*cos(theta1+theta2) ...
     + m4*g*lc3*cos(theta1+theta2+theta3);

G3 = m4*g*lc3*cos(theta1+theta2+theta3);

G = [G1; G2; G3];
end
