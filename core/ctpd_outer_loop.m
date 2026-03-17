function tau_des = ctpd_outer_loop(q, dq, q_des, dq_des, ddq_des, P, Phi_k)
% CTPD_OUTER_LOOP  Computed-torque PD outer loop (Eqs. 21-22 in paper).
%
%   Computes the desired joint torque using full nonlinear dynamics
%   compensation plus PD error feedback.
%
%   tau_des = M(q)*ddq_cmd + C(q,dq)*dq + G(q) + D*dq
%   ddq_cmd = ddq_des + Kp*(q_des - q) + Kd*(dq_des - dq)
%
%   Inputs:
%     q       - n x 1  current joint angles [rad]
%     dq      - n x 1  current joint velocities [rad/s]
%     q_des   - n x 1  desired joint angles [rad]
%     dq_des  - n x 1  desired joint velocities [rad/s]
%     ddq_des - n x 1  desired joint accelerations [rad/s^2]
%     P       - system parameter struct
%     Phi_k   - nw x nw deterioration matrix (unused in outer loop,
%               included for interface consistency)
%
%   Output:
%     tau_des - n x 1  desired joint torque [Nm]

%% PD commanded acceleration
e   = q_des  - q;    % position error
de  = dq_des - dq;   % velocity error

ddq_cmd = ddq_des + P.Kp * e + P.Kd * de;

%% Nonlinear dynamics terms
M_q = sts_inertia_public(q, P);
C_q = sts_coriolis_public(q, dq, P);
G_q = sts_gravity_public(q, P);

%% Desired torque
tau_des = M_q * ddq_cmd + (C_q + P.D_damp) * dq + G_q;

%% Saturate to motor torque limits
tau_lim = P.tau_max_motor;
tau_des = max(-tau_lim, min(tau_lim, tau_des));

end


% ─── helper wrappers (calls shared functions from sts_dynamics.m) ────────

function M = sts_inertia_public(q, P)
% Thin wrapper so ctpd can call the same inertia model.
% (In practice, consolidate into a single shared utility file.)
m  = P.m_seg;
l  = P.l_seg;
lc = l / 2;
I  = (1/12) * m .* l.^2;

theta2 = q(2);
theta3 = q(3);

m2 = m(2); m3 = m(3); m4 = m(4);
l2 = l(2); l3 = l(3);
lc3 = lc(3);
I2 = I(2); I3 = I(3); I4 = I(4);

M11 = (m2+m3+m4)*l2^2 + (m3+m4)*l3^2 + I2+I3+I4 ...
    + 2*(m3+m4)*l2*l3*cos(theta2) ...
    + 2*m4*l2*lc3*cos(theta2) + 2*m4*l3*lc3*cos(theta3);
M12 = (m3+m4)*l3^2 + I3+I4 ...
    + (m3+m4)*l2*l3*cos(theta2) + m4*l3*lc3*cos(theta3);
M13 = m4*lc3^2 + I4;
M22 = (m3+m4)*l3^2 + I3+I4 + 2*m4*l3*lc3*cos(theta3);
M23 = m4*lc3^2 + I4 + m4*l3*lc3*cos(theta3);
M33 = m4*lc3^2 + I4;

M = [M11 M12 M13; M12 M22 M23; M13 M23 M33];
M = 0.5*(M+M') + 1e-6*eye(3);
end

function C = sts_coriolis_public(q, dq, P)
l = P.l_seg; m = P.m_seg;
theta2 = q(2); theta3 = q(3);
m3 = m(3); m4 = m(4);
l2 = l(2); l3 = l(3);
h2 = -(m3+m4)*l2*l3*sin(theta2);
h3 = -m4*l2*l3*sin(theta3);
C = [h2*(2*dq(2)+dq(3))+h3*dq(3), h2*dq(1), h3*dq(1);
     -h2*dq(1), 0, h3*(dq(1)+dq(2));
     -(h3*(2*dq(1)+dq(2))), -h3*dq(2), 0];
end

function G = sts_gravity_public(q, P)
m = P.m_seg; l = P.l_seg; lc = l/2; g = 9.81;
m2=m(2); m3=m(3); m4=m(4);
l2=l(2); l3=l(3); lc3=lc(3);
theta1=q(1); theta2=q(2); theta3=q(3);
G1 = (m2+m3+m4)*g*l2*cos(theta1) + (m3+m4)*g*l3*cos(theta1+theta2) ...
     + m4*g*lc3*cos(theta1+theta2+theta3);
G2 = (m3+m4)*g*l3*cos(theta1+theta2) + m4*g*lc3*cos(theta1+theta2+theta3);
G3 = m4*g*lc3*cos(theta1+theta2+theta3);
G  = [G1; G2; G3];
end
