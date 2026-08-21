function [L, g2] = Lemma1(System)
% Lemma 1: H-infinity Observer Gain Design via LMI
%
% State estimation error dynamics (single agent):
%   ebar(k+1) = (A - L*C)*ebar(k) + E*w(k) - L*F*v(k)
%
% H-infinity criterion:
%   sum_{k} ||ebar(k)||^2 < gamma^2 * sum_{k} (||w(k)||^2 + ||v(k)||^2)
%
% LMI formulation (Schur complement + substitution Y = P*L):
%   [ -P+I,       0,        0,      A'P - C'Y' ]
%   [   0,    -g2*I,        0,          E'P    ] < 0,  P > 0
%   [   0,        0,    -g2*I,         -F'Y'   ]
%   [ PA-YC,    PE,      -YF,            -P    ]
%
% Observer gain recovered as: L = P^{-1} * Y
% H-infinity performance index: gamma = sqrt(g2)

A  = System.A;
C  = System.C;
E  = System.E;
F  = System.F;

n  = size(A,1);
q  = size(C,1);
nw = size(E,2);
nv = size(F,2);

eps = 1e-10;

%% Variable Setup
P  = sdpvar(n, n, 'symmetric');
Y  = sdpvar(n, q, 'full');
g2 = sdpvar(1);

%% LMI Constraints
LMI = [];
LMI = [LMI, P >= eps * eye(n)];
LMI = [LMI, g2 >= 0];

M11 = -P + eye(n);
M14 = A'*P - C'*Y';
M22 = -g2 * eye(nw);
M24 = E'*P;
M33 = -g2 * eye(nv);
M34 = -F'*Y';
M41 = P*A - Y*C;
M42 = P*E;
M43 = -Y*F;
M44 = -P;

M = [M11,           zeros(n, nw),  zeros(n, nv),  M14;
     zeros(nw, n),  M22,           zeros(nw, nv), M24;
     zeros(nv, n),  zeros(nv, nw), M33,           M34;
     M41,           M42,           M43,           M44];

LMI = [LMI, M <= -eps * eye(size(M, 1))];

%% Solve: minimize gamma^2
Objective = g2;
options   = sdpsettings('solver', 'mosek', 'verbose', 0);
sol       = optimize(LMI, Objective, options);

if sol.problem == 0
    P_val = value(P);
    Y_val = value(Y);
    g2    = value(g2);
    L     = P_val \ Y_val;

    % Check: eigenvalues of (A - L*C) must be inside the unit circle
    eig_obs = eig(A - L*C);
    if all(abs(eig_obs) < 1)
        fprintf('Lemma 1 solved successfully. gamma = %.6f\n', sqrt(g2));
        fprintf('Observer (A-LC) is stable: max|eig| = %.6f\n', max(abs(eig_obs)));
    else
        warning('Lemma 1: LMI solved but (A-LC) has eigenvalues outside the unit circle.');
        fprintf('Eigenvalues of (A-LC): \n'); disp(eig_obs);
    end
else
    error('Lemma 1: LMI infeasible or solver error. %s', sol.info);
end

end
