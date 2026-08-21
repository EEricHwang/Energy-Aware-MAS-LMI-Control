function [K, g1, P0_val] = Theorem1(System)
% Theorem 1: H-infinity Consensus Controller Gain Design via LMI
%
% Consensus error dynamics (estimation error neglected; gamma from Lemma 1 is small):
%   δ̃_i(k+1) = (A - ρ_s*λ_i*B*K)*δ̃_i(k) + E*w̃_i(k)
%
% where:
%   λ_i  : i-th non-zero eigenvalue of L_graph  (i = 2,...,N)
%   ρ_s  : SOC-weight vertex  ρ ∈ {w_under, w_upper}
%   w_under = exp(-2*alpha/soc_min),  w_upper = exp(-2*alpha/soc_max)
%
% SOC-dependent weight: w_ij = exp(-alpha*(1/soc_i + 1/soc_j))
%   w_under = exp(-2*alpha/soc_min),  w_upper = exp(-2*alpha/soc_max)
%
% Polytopic (per-mode) LMI — substitution Y = K*P0, P0 = P^{-1}:
%   For each vertex s and mode i:
%   [ -P0,    0,         (A*P0 - ρ_s*λ_i*B*Y)',   P0 ]
%   [  0,    -g1*I,       E',                       0  ] < 0
%   [ A*P0-ρ_s*λ_i*B*Y,  E,   -P0,                 0  ]
%   [  P0,    0,           0,                      -I  ]
%
% Controller gain recovered as: K = Y * P0^{-1}
% H-infinity performance index: gamma = sqrt(g1)
%
% Note: System must include fields: A, B, E, L, soc_min, soc_max, alpha

A      = System.A;
B      = System.B;
E      = System.E;
L_topo = System.L;

soc_min = System.soc_min;
soc_max = System.soc_max;
alpha   = System.alpha;

n  = size(A, 1);
m  = size(B, 2);
nw = size(E, 2);

w_under = exp(-2*alpha / soc_min);
w_upper = exp(-2*alpha / soc_max);
rho     = [w_under, w_upper];

% Laplacian eigendecomposition (modal decomposition)
[~, Lam]         = eig(L_topo);
lambda_all       = real(diag(Lam));
lambda_all       = sort(lambda_all);
lambda_c         = lambda_all(2:end);   % non-zero consensus modes (N-1)

fprintf('=== Graph Eigenvalues ===\n');
fprintf('  λ_c = '); fprintf('%.4f  ', lambda_c); fprintf('\n');
fprintf('  w_under = %.6f,  w_upper = %.6f\n\n', w_under, w_upper);

eps    = 1e-10;
sz_lmi = n + nw + n + n;   % 4+4+4+4 = 16

%% Variable Setup
P0 = sdpvar(n, n, 'symmetric');    % P0 = P^{-1}  (n×n)
Y  = sdpvar(m, n, 'full');         % Y  = K*P0    (m×n)
g1 = sdpvar(1);

%% LMI Constraints
LMI = [];
LMI = [LMI, P0 >= eps * eye(n)];
LMI = [LMI, g1 >= eps];

% Polytopic vertices × consensus modes
for s = 1:2
    for i = 1:length(lambda_c)
        rho_s = rho(s);
        lam_i = lambda_c(i);

        AsP0 = A*P0 + rho_s*lam_i*(B*Y);   % (A - ρ_s*λ_i*B*K)*P0

        Mi = [-P0,              zeros(n, nw),       AsP0',              P0;
               zeros(nw, n),   -g1*eye(nw),         E',                 zeros(nw, n);
               AsP0,            E,                  -P0,                zeros(n);
               P0,              zeros(n, nw),        zeros(n),          -eye(n)];

        LMI = [LMI, Mi <= -eps * eye(sz_lmi)];
    end
end

%% Solve: minimize gamma^2
Objective = g1;
options   = sdpsettings('solver', 'mosek', 'verbose', 0, 'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-8);
sol       = optimize(LMI, Objective, options);

if sol.problem == 0
    P0_val = value(P0);
    Y_val  = value(Y);
    g1     = value(g1);
    K      = Y_val / P0_val;    % K = Y * P0^{-1}

    fprintf('Theorem 1 solved successfully. gamma = %.6f\n', sqrt(g1));

    % Check: eigenvalues of (A + B*K) inside unit circle
    eig_ctrl = eig(A + B*K);
    if all(abs(eig_ctrl) < 1)
        fprintf('Controller (A+BK) is stable: max|eig| = %.6f\n', max(abs(eig_ctrl)));
    else
        warning('Theorem 1: LMI solved but (A+BK) has eigenvalues outside the unit circle.');
        fprintf('Eigenvalues of (A+BK): \n'); disp(eig_ctrl);
    end

    % Verify all closed-loop modes across vertices
    fprintf('\n--- Closed-loop verification (all modes & vertices) ---\n');
    all_ok = true;
    for s = 1:2
        for i = 1:length(lambda_c)
            A_si = A + rho(s)*lambda_c(i)*B*K;
            ev   = max(abs(eig(A_si)));
            ok   = ev < 1;
            all_ok = all_ok & ok;
            fprintf('  Vertex %d, Mode %d (rho=%.4f, lambda=%.4f): max|eig|=%.6f  %s\n', s, i, rho(s), lambda_c(i), ev, status_str(ok));
        end
    end
    if all_ok
        fprintf('  >> All modes STABLE\n');
    else
        warning('Theorem 1: Some closed-loop modes are UNSTABLE.');
    end
else
    error('Theorem 1: LMI infeasible or solver error. %s', sol.info);
end

end

%% Helper
function s = status_str(ok)
    if ok; s = 'STABLE'; else; s = 'UNSTABLE'; end
end
