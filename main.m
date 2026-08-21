clear all
close all
clc

%% 1. MAS Setup (XY Double integrator)
dT       = 0.05;                                        % [-] sampling time
N        = 5;                                           % [-] number of agents
choice   = 1;                                           % [-] network topology number

A        = [1 dT 0 0; 0 1 0 0; 0 0 1 dT; 0 0 0 1];      % [-] system matrix
B        = [dT^2/2 0; dT 0; 0 dT^2/2; 0 dT];            % [-] input matrix
C        = [1 0 0 0; 0 0 1 0];                          % [-] output matrix
E        = 10e-5*eye(4);                                % [-] process noise matrix
F        = 10e-5*eye(2);                                % [-] measurement noise matrix

n        = size(A,1); 
m        = size(B,2);
q        = size(C,1);
tilde_w  = size(E,1);
tilde_v  = size(F,1);
M        = eye(N) - (1/N)*ones(N,1)*ones(N,1)';

soc_min  = 0.2;
soc_max  = 0.9;
alpha    = 0.5;

w_under  = exp(-(2*alpha)/soc_min);
w_upper  = exp(-(2*alpha)/soc_max);
w_star   = w_upper - w_under;

switch choice
    case 1
        L = [ 3 -1 -1 -1  0;   % Agent 1
             -1  2 -1  0  0;   % Agent 2
             -1 -1  4 -1 -1;   % Agent 3
             -1  0 -1  3 -1;   % Agent 4
              0  0 -1 -1  2];  % Agent 5
    case 2
        L = [ 2 -1 -1; 
             -1  2 -1;
             -1 -1  2];
end

System.A        = A;
System.B        = B;
System.C        = C;
System.E        = E;
System.F        = F;
System.L        = L;
System.M        = M;
System.soc_min  = soc_min;
System.soc_max  = soc_max;
System.alpha    = alpha;
System.Ts       = dT;

% Battery model parameters (used in Theorem 2 & Simulation)
System.eta_ctrl = 1e-9;    % SOC drain per unit ||u_i||^2
System.eta_prop = 5e-5;    % SOC drain per unit ||v_i||^2 (propulsion)

% Theorem 2 design parameters
System.e_max    = 2.0;     % initial consensus error bound: ||e(0)|| <= e_max
System.v_max    = 0.5;     % maximum agent speed bound: ||v_i|| <= v_max

%% Lemma 1: Observer gain
[L_obs, g2] = Lemma1(System);

%% Theorem 1: Controller gain
[K, g1, P0_val] = Theorem1(System);

%% Theorem 2: Battery safety guarantee
%         1;    2;   3;   4;   5
soc0 = [0.9; 0.85; 0.7; 0.6; 0.5];    % initial SOC (same as Simulation)

N = size(L, 1);

% SOC-dependent weight matrix
W = zeros(N, N);
for i = 1:N
    for j = 1:N
        if i ~= j && L(i,j) ~= 0
            W(i,j) = exp(-alpha * (1/soc0(i) + 1/soc0(j)));
        end
    end
end

[u_max, C_rate_max, T_safe] = Theorem2(System, K, P0_val, soc0);

%% Simulation
Simulation(System, L_obs, K, u_max, C_rate_max, T_safe);

%% Export gains to YAML (for ROS2 & Gazebo Implementation)
export_gains(System, L_obs, K, g1, g2, dT);
