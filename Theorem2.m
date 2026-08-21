function [u_max, C_rate_max, T_safe] = Theorem2(System, K, P0_val, soc0)
% Theorem 2: Battery Safety Guarantee via Control Input Bounding
%
% Under the controller K designed in Theorem 1, this function verifies that:
%   (1) ||u_i(k)||^2  <= u_max^2      for all k >= 0, for all i
%   (2) C-rate_i(k)   <= C_rate_max   for all k >= 0, for all i
%   (3) soc_i(k)      >= soc_min      for all k <= T_safe, for all i
%
% Key idea:
%   From Theorem 1, the consensus error trajectory stays inside the
%   Lyapunov level set:
%       e^T(k) (I_N x P) e(k) <= c_bar,  for all k >= 0
%
%   where c_bar = lambda_max(P) * e_max^2  (worst-case initial condition).
%
%   The LMI condition (Schur complement):
%       [ u_max^2/c_bar * I_m,   Psi_i_max           ]
%       [ Psi_i_max^T,           I_N (x) P           ] >= 0
%
%   guarantees ||u_i(k)||^2 <= u_max^2 for all e in the level set.
%
%   Here, Psi_i_max = K * [w_upper * L (x) I_n]_i  (worst-case weight row)
%
%   Battery safety then follows from the SOC dynamics:
%       soc_i(k+1) = soc_i(k) - eta_ctrl*||u_i(k)||^2 - eta_prop*||v_i(k)||^2
%                 >= soc_i(k) - delta_max
%
% Inputs:
%   System   - system struct with fields:
%                A, B, L, alpha, soc_min, soc_max, Ts
%                eta_ctrl  [SOC/control_energy]
%                eta_prop  [SOC/speed_energy]
%                e_max     initial consensus error bound: ||e(0)|| <= e_max
%                v_max     maximum agent speed bound:    ||v_i||  <= v_max
%   K        - controller gain (m x n) from Theorem 1
%   P0_val   - P0 = P^{-1} (n x n) from Theorem 1
%   soc0     - initial SOC vector (N x 1)
%
% Outputs:
%   u_max      - guaranteed control input bound  sqrt(u_max^2)
%   C_rate_max - guaranteed C-rate upper bound   [SOC/step]
%   T_safe     - safe operation horizon          [steps]

%% ── Unpack ──────────────────────────────────────────────────────────────
L_topo   = System.L;
alpha    = System.alpha;
soc_min  = System.soc_min;
soc_max  = System.soc_max;
eta_ctrl = System.eta_ctrl;
eta_prop = System.eta_prop;
Ts       = System.Ts;
e_max    = System.e_max;
v_max    = System.v_max;

N        = size(L_topo, 1);
n        = size(K, 2);
m        = size(K, 1);

%% ── Derived quantities ──────────────────────────────────────────────────
P        = inv(P0_val);                          % Lyapunov matrix  (n x n), 역행렬 씌우기
P_big    = kron(eye(N), P);                      % I_N (x) P        (nN x nN)
w_upper  = exp(-2 * alpha / soc_max);            % worst-case SOC weight

% c_bar: upper bound on e^T(0) (I_N x P) e(0)
%   e^T(I_N x P) e <= lambda_max(P) * ||e||^2 <= lambda_max(P) * e_max^2
c_bar    = max(real(eig(P))) * e_max^2;

fprintf('=== Theorem 2: Battery Safety Analysis ===\n');
fprintf('  P eigenvalue range : [%.4f, %.4f]\n', min(real(eig(P))), max(real(eig(P))));
fprintf('  c_bar (level set)  : %.6f\n', c_bar);
fprintf('  w_upper            : %.6f\n', w_upper);
fprintf('  e_max              : %.4f\n', e_max);
fprintf('  v_max              : %.4f\n\n', v_max);

%% ── LMI: find minimum u_max^2 satisfying Schur condition ───────────────
%
%   Schur complement condition for each agent i:
%       [ u2/c_bar * I_m,   Psi_i_max ] >= 0
%       [ Psi_i_max^T,      P_big     ]
%
%   Equivalent to:
%       u2 >= c_bar * lambda_max( Psi_i_max * (I_N x P0) * Psi_i_max^T )
%
eps_lmi  = 1e-10;
sz_schur = m + N * n;

u2       = sdpvar(1);
LMI      = [];
LMI      = [LMI, u2 >= eps_lmi];

for i = 1:N
    % i-th block row of (w_upper * L (x) I_n):  size n x nN
    L_row_i  = L_topo(i, :);
    Psi_i    = K * kron(w_upper * L_row_i, eye(n));     % m x nN

    M_schur  = [ u2 / c_bar * eye(m),   Psi_i  ;
                 Psi_i',                P_big  ];

    LMI = [LMI, M_schur >= eps_lmi * eye(sz_schur)];
end

options = sdpsettings('solver', 'mosek', 'verbose', 0, 'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-8);
sol     = optimize(LMI, u2, options);

if sol.problem ~= 0
    error('Theorem 2: LMI infeasible or solver error. %s', sol.info);
end

u_max2  = value(u2);
u_max   = sqrt(u_max2);

%% ── C-rate bound ─────────────────────────────────────────────────────────
%   C-rate_i(k) = eta_ctrl * ||u_i(k)||^2 <= eta_ctrl * u_max^2
C_rate_max = eta_ctrl * u_max2;     % [SOC/step]

%% ── SOC safety horizon ───────────────────────────────────────────────────
%   soc_i(k) >= soc_i(0) - k * delta_max >= soc_min
%   => k <= (soc_i(0) - soc_min) / delta_max  =: T_safe
delta_max = eta_ctrl * u_max2 + eta_prop * v_max^2;
T_safe    = floor((min(soc0) - soc_min) / delta_max);

%% ── Print results ────────────────────────────────────────────────────────
fprintf('--- Theorem 2 Results ---\n');
fprintf('  u_max          = %.6f\n',     u_max);
fprintf('  C_rate_max     = %.4e  [SOC/step]\n', C_rate_max);
fprintf('  delta_max      = %.4e  [SOC/step]\n', delta_max);
fprintf('  T_safe         = %d steps  (%.1f s)\n', T_safe, T_safe * Ts);
fprintf('  min soc0       = %.2f,  soc_min = %.2f\n', min(soc0), soc_min);

if T_safe > 0
    fprintf('  >> Battery safety GUARANTEED for %d steps (%.1f s)\n\n', T_safe, T_safe * Ts);
else
    warning('Theorem 2: T_safe <= 0. Increase e_max margin or check soc0.');
end

%% ── Cross-check via direct eigenvalue computation ────────────────────────
%   u_max^2 = max_i { c_bar * lambda_max( Psi_i * (I_N x P0) * Psi_i^T ) }
P0_big   = kron(eye(N), P0_val);
u2_check = 0;
for i = 1:N
    L_row_i  = L_topo(i, :);
    Psi_i    = K * kron(w_upper * L_row_i, eye(n));
    M_eig    = Psi_i * P0_big * Psi_i';
    u2_check = max(u2_check, c_bar * max(real(eig(M_eig))));
end
fprintf('  Cross-check (eigenvalue): u_max = %.6f  |  LMI: %.6f\n', sqrt(u2_check), u_max);
if abs(sqrt(u2_check) - u_max) < 1e-4
    fprintf('  >> LMI and eigenvalue results CONSISTENT\n\n');
else
    warning('Theorem 2: LMI and eigenvalue results differ. Check numerics.');
end

end
