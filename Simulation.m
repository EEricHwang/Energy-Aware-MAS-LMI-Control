function Simulation(System, L_obs, K, u_max, C_rate_max, T_safe)
% Simulation: Energy-Aware Pentagon Formation Control
%
% 5 agents maintain a regular pentagon formation while drifting
% diagonally in the x+ y+ direction.
%
% Formation control law:
%   u_i(k) = K * sum_{j in N_i} w_ij * (zhat_i(k) - zhat_j(k))
%   where zhat_i = xhat_i - delta_i  (shifted state = formation error state)
%
% SOC drain model (two sources):
%   soc_i(k+1) = soc_i(k) - eta_ctrl * ||u_i(k)||^2        (control effort)
%                         - eta_prop * (vx_i^2 + vy_i^2)   (propulsion energy)
%
% => SOC keeps decreasing even after formation is achieved (because agents
%    are still moving), but the proposed method is less aggressive.
%
% Inputs (Theorem 2 outputs — battery safety bounds):
%   u_max      : guaranteed control input bound from Theorem 2
%   C_rate_max : guaranteed C-rate upper bound  [SOC/step]
%   T_safe     : safe operation horizon         [steps]
%
% Causal order per step k:
%   1. Compute u(k) from zhat(k)
%   2. Measure y(k) from x(k)   [before plant update]
%   3. Update xhat(k+1) = A*xhat + B*u + L_obs*(y - C*xhat)
%   4. Update x(k+1) = A*x + B*u + E*w

Case     = 1;

%% Unpack system
A        = System.A;
B        = System.B;
C        = System.C;
E        = System.E;
F        = System.F;
L_topo   = System.L;
alpha    = System.alpha;
soc_min  = System.soc_min;
soc_max  = System.soc_max;
eta_ctrl = System.eta_ctrl;
eta_prop = System.eta_prop;
N        = size(L_topo, 1);    % 5
n        = size(A, 1);         % 4: [px; vx; py; vy]
q        = size(C, 1);         % 2

%% Pentagon formation offsets (position only; velocity offset = 0)
R = 2.0;                                           % formation radius [m]
theta = 2*pi*(0:N-1)'/N + pi/2;                    % top vertex first
delta = zeros(n, N);
for i = 1:N
    delta(:,i) = [R*cos(theta(i)); 0; R*sin(theta(i)); 0];
end

% Diagonal reference velocity (x+ y+ direction, 45 deg)
v_ref  = 0.3;                               % formation speed [m/s]
vx_ref = v_ref / sqrt(2);
vy_ref = v_ref / sqrt(2);

% Simulation parameters
T_sim  = 20000;
dT     = System.Ts;
t      = (0:T_sim-1)*dT;

% Battery model
soc0 = [0.9; 0.85; 0.7; 0.6; 0.5];          % heterogeneous initial SOC

%% Initial conditions: agents near their pentagon positions, moving diagonally
switch Case
    case 1
        % rng(42);
        % x0 = delta;
        % for i = 1:N
        %     x0(1,i) = x0(1,i) + 0.5*randn;     % perturb px
        %     x0(2,i) = vx_ref + 0.05*randn;     % vx near reference
        %     x0(3,i) = x0(3,i) + 0.5*randn;     % perturb py
        %     x0(4,i) = vy_ref + 0.05*randn;     % vy near reference
        % end
        % xhat0 = x0 + 0.1*randn(n, N);

        x0(1,1) = -0.2691;
        x0(2,1) = 0.2555;
        x0(3,1) = 2.4880;
        x0(4,1) = 0.2290;

        x0(1,2) = -2.4002;
        x0(2,2) = 0.1860;
        x0(3,2) = -0.0307;
        x0(4,2) = 0.2580;

        x0(1,3) = -1.0873;
        x0(2,3) = 0.2499;
        x0(3,3) = -1.9138;
        x0(4,3) = 0.3044;

        x0(1,4) = 2.0840;
        x0(2,4) = 0.2059;
        x0(3,4) = -2.1733;
        x0(4,4) = 0.1781;

        x0(1,5) = 1.9092;
        x0(2,5) = 0.2092;
        x0(3,5) = 0.2875;
        x0(4,5) = 0.2274;

        %%
        xhat0(1,1) = -0.3100;
        xhat0(2,1) = 0.1273;
        xhat0(3,1) = 2.4595;
        xhat0(4,1) = 0.2225;

        xhat0(1,2) = -2.3001;
        xhat0(2,2) = 0.1058;
        xhat0(3,2) = -0.0215;
        xhat0(4,2) = 0.2854;

        xhat0(1,3) = -1.3006;
        xhat0(2,3) = 0.2903;
        xhat0(3,3) = -2.0096;
        xhat0(4,3) = 0.2717;

        xhat0(1,4) = 2.3319;
        xhat0(2,4) = 0.3598;
        xhat0(3,4) = -2.0906;
        xhat0(4,4) = 0.1621;

        xhat0(1,5) = 1.8234;
        xhat0(2,5) = 0.2490;
        xhat0(3,5) = 0.2650;
        xhat0(4,5) = 0.0960;
    case 2
        % state: [px; vx; py; vy]
        vx_ref = v_ref / sqrt(2);   % ≈ 0.2121 m/s
        vy_ref = v_ref / sqrt(2);   % ≈ 0.2121 m/s
        
        x0 = [ 0.3,  -0.3,   0.2,  -0.2,   0.1;          % px
               vx_ref, vx_ref, vx_ref, vx_ref, vx_ref;   % vx
               0.2,   0.3,  -0.3,  -0.1,   0.1;          % py
               vy_ref, vy_ref, vy_ref, vy_ref, vy_ref];  % vy
        
        xhat0 = x0;   
end

% Shared noise realizations
rng(42);
sigma_w = 10e-5;
sigma_v = 10e-5;
W_noise = sigma_w * randn(n, N, T_sim);
V_noise = sigma_v * randn(q, N, T_sim);

% Run both methods
res1 = run_sim(A, B, C, E, F, L_topo, K, L_obs, delta, ...
    x0, xhat0, soc0, W_noise, V_noise, ...
    T_sim, alpha, eta_ctrl, eta_prop, soc_min, soc_max, 'proposed');

res2 = run_sim(A, B, C, E, F, L_topo, K, L_obs, delta, ...
    x0, xhat0, soc0, W_noise, V_noise, ...
    T_sim, alpha, eta_ctrl, eta_prop, soc_min, soc_max, 'Baseline');

%% 논문에 들어갈 그림 선별 
lw   = 2;
cols = lines(N);
leg_agents = arrayfun(@(i) sprintf('Agent %d (soc_0=%.1f)', i, soc0(i)), 1:N, 'UniformOutput', false);

T_show = 500;   
cols = lines(N);

figure(5)
set(gcf,'Color','w','Position',[200 200 1000 800])
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- (a) XY Trajectory ---
nexttile
h_prop = gobjects(N,1);
for i = 1:N
    h_prop(i) = plot(res1.px(i, 1:T_show), res1.py(i, 1:T_show), '-',  'Color', cols(i,:), 'LineWidth', 1.5); hold on;
    plot(res2.px(i, 1:T_show), res2.py(i, 1:T_show), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    plot(res1.px(i,1), res1.py(i,1), 'o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end

% Pentagon formation snapshots
pent_ord    = [1 2 3 4 5 1];
snap_ks     = [1, round(T_show/2), T_show];
alphas      = [0.30, 0.60, 1.00];
lws_snap    = [0.8,  1.2,  2.0];
snap_labels = {'Formation (t=0)', ...
               sprintf('Formation (t=%.0fs)', t(snap_ks(2))), ...
               sprintf('Formation (t=%.0fs)', t(snap_ks(3)))};
h_snap = gobjects(3,1);
for s = 1:3
    k = snap_ks(s);
    h_snap(s) = plot(res1.px(pent_ord,k), res1.py(pent_ord,k), 'k-', ...
        'LineWidth', lws_snap(s), 'Color', [0 0 0 alphas(s)]);
end
for i = 1:N
    plot(res1.px(i,T_show), res1.py(i,T_show), 's', ...
        'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
end
h_conv = plot(nan, nan, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
leg_entries = [h_prop(:); h_conv; h_snap(1); h_snap(3)];
leg_labels  = [arrayfun(@(i) sprintf('Agent %d (Proposed)', i), 1:N, 'UniformOutput', false), ...
               {'Baseline (all agents)'}, snap_labels(1), snap_labels(3)];
legend(leg_entries, leg_labels, 'Location', 'northwest', 'FontSize', 8);

xlabel('X Position [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position [m]', 'FontSize', 12, 'FontWeight', 'bold');
title('XY Trajectory for [0 25]s: Proposed vs. Baseline', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-3 8.5]); ylim([-3 8.5]);
axis equal; grid on;

% --- (b) SOC comparison ---
nexttile
semilogy(t, res1.form_err, 'b-',  'LineWidth', lw); hold on;
semilogy(t, res2.form_err, 'r--', 'LineWidth', lw);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\Sigma_i ||\rho_i - \bar{\rho}||_2$', 'Interpreter','latex','FontSize', 18, 'FontWeight', 'bold');
title('Formation Error','FontSize', 12,'FontWeight','bold');
legend('Proposed (Energy-Aware)', 'Baseline', 'Location','northeast','FontSize',10);
grid on;

% --- (c) Control input for Proposed + u_max ---
nexttile
control_prop = gobjects(N,1);
for i = 1:N
    control_prop(i) = plot(t, res1.u_each(i,:), '-',  'Color', cols(i,:), 'LineWidth', 2.5); hold on;
end

annotation('textbox', [0.08 0.315 0.2 0.1], ...
    'String', 'No Safety Violations', ...
    'FontSize', 16, ...
    'Color', 'b', ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none');

yline(u_max, 'k-.', 'LineWidth', 2.0);
text(15, u_max*0.89, '$u_{\max}$ (Theorem 2)', 'Interpreter', 'latex', 'FontSize', 14, 'VerticalAlignment', 'bottom');
leg_entries = [control_prop(:)];
leg_labels  = [arrayfun(@(i) sprintf('Agent %d (Proposed)', i), 1:N, 'UniformOutput', false)];
legend(leg_entries, leg_labels, 'Location', 'southeast', 'FontSize', 10);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\|u_i(k)\|_2$', 'Interpreter','latex','FontSize', 18, 'FontWeight', 'bold');
title('Control Input: Proposed', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T_show*dT])
grid on;

% --- (d) Control input for Baseline + u_max ---
nexttile
control_base = gobjects(N,1);
for i = 1:N
    control_base(i) = plot(t, res2.u_each(i,:), '--', 'Color', cols(i,:), 'LineWidth', 2.5); hold on;
end

annotation('arrow', [0.59 0.55], [0.22 0.26], 'Color', 'r', 'LineWidth', 2);
annotation('textbox', [0.59 0.135 0.2 0.1], ...
    'String', 'Safety Violations', ...
    'FontSize', 16, ...
    'Color', 'r', ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none');

yline(u_max, 'k-.', 'LineWidth', 2.0);
text(15, u_max*1.04, '$u_{\max}$ (Theorem 2)', 'Interpreter','latex','FontSize',14,'VerticalAlignment','bottom');
leg_entries = [control_base(:)];
leg_labels  = [arrayfun(@(i) sprintf('Agent %d (Baseline)', i), 1:N, 'UniformOutput', false)];
legend(leg_entries, leg_labels, 'Location', 'northeast', 'FontSize', 10);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$\|u_i(k)\|_2$', 'Interpreter','latex','FontSize', 18, 'FontWeight', 'bold');
title('Control Input: Baseline', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T_show*dT])
grid on;

%%
figure(7)
set(gcf,'Color','w','Position',[200 200 1000 400])
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- (a) Cumulative energy consumption ---
% nexttile
% plot(t, res1.cum_energy, 'b-',  'LineWidth', 2.0); hold on;
% plot(t, res2.cum_energy, 'r--', 'LineWidth', 2.0);
% xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('$\sum_i \|u_i(k)\|^2$ (cumulative)', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold');
% title('Cumulative Energy Consumption', 'FontSize', 12, 'FontWeight', 'bold');
% legend('Proposed', 'Baseline', 'Location', 'northwest', 'FontSize', 10);
% grid on;

nexttile
for i = 1:N
    plot(t, res1.soc(i,:), '-',  'Color', cols(i,:), 'LineWidth', 2.5); hold on;
    plot(t, res2.soc(i,:), '--', 'Color', cols(i,:), 'LineWidth', 2.5);
end
t_safe_sec = T_safe * dT;          
xline(t_safe_sec, 'm-.', 'LineWidth', 2.0);   
text(490, 0.11, '$T_{\mathrm{safe}} = 843.4$ s', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'm', 'VerticalAlignment', 'bottom');
yline(soc_min, 'k--', 'LineWidth', 1.5);
yline(soc_max, 'k--', 'LineWidth', 1.5);
text(t(end)*0.02, soc_max+0.02, '$\mathrm{soc}_{\max} = 0.9$', 'Interpreter', 'latex', 'FontSize', 14, 'VerticalAlignment', 'bottom');
text(t(end)*0.02, soc_min+0.02, '$\mathrm{soc}_{\min} = 0.2$', 'Interpreter', 'latex', 'FontSize', 14, 'VerticalAlignment', 'bottom');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('SOC', 'FontSize', 12, 'FontWeight', 'bold');
title('SOC Trajectories', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0.1 1]); grid on;

% --- (b) SOC-adaptive weight evolution (Proposed only) ---
% Extract undirected edges from L_topo
edges = [];
for i = 1:N
    for j = i+1:N
        if L_topo(i,j) ~= 0
            edges = [edges; i j];
        end
    end
end
n_edges = size(edges, 1);

% Compute w_ij(k) from stored SOC trajectories (Proposed)
w_traj = zeros(n_edges, T_sim);
for k = 1:T_sim
    soc_k = max(res1.soc(:,k), soc_min);
    for e = 1:n_edges
        ii = edges(e,1);
        jj = edges(e,2);
        w_traj(e,k) = exp(-alpha * (1/soc_k(ii) + 1/soc_k(jj)));
    end
end

% Reference bounds
w_upper_ref = exp(-2*alpha / soc_max);   % max weight (both agents at soc_max)
w_lower_ref = exp(-2*alpha / soc_min);   % min weight (both agents at soc_min)

nexttile
edge_cols = lines(n_edges);
h_w = gobjects(n_edges, 1);
for e = 1:n_edges
    h_w(e) = plot(t, w_traj(e,:), '-', 'Color', edge_cols(e,:), 'LineWidth', 2.5); hold on;
end
yline(w_upper_ref, 'k--', 'LineWidth', 1.5);
yline(w_lower_ref, 'k--',  'LineWidth', 1.5);
text(t(end)*0.7, w_upper_ref*1.02, '$\bar{w} = 0.329$', 'Interpreter','latex','FontSize',13,'VerticalAlignment','bottom');
text(t(end)*0.7, w_lower_ref*2.1, '$\underline{w} = 6.74e^{-2}$', 'Interpreter','latex','FontSize',13,'VerticalAlignment','bottom');
leg_labels_w = arrayfun(@(e) sprintf('$w_{%d%d}$', edges(e,1), edges(e,2)), 1:n_edges, 'UniformOutput', false);
legend(h_w, leg_labels_w, 'Interpreter','latex', 'Location','southwest', 'FontSize', 11);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('$w_{ij}(k)$', 'Interpreter','latex', 'FontSize', 14, 'FontWeight', 'bold');
title('SOC-Dependent Weight Evolution', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 0.38])
grid on;

end

%% ======================================================
%  Battery safety check helper
%% ======================================================
function s = check_str(ok)
    if ok; s = 'YES ✓'; else; s = 'NO  ✗'; end
end

%% ======================================================
%  Pentagon snapshot helper
%% ======================================================
function draw_pentagon(px, py, snap_idx, cols, N, mk_size)
for s = 1:length(snap_idx)
    k = snap_idx(s);
    alpha_val = 0.3 + 0.7*(s-1)/(length(snap_idx)-1);  % darker = later
    pts_x = [px(:,k); px(1,k)];
    pts_y = [py(:,k); py(1,k)];
    plot(pts_x, pts_y, 'k-', 'LineWidth', 1.0, 'Color', [0 0 0 alpha_val]); hold on;
    for i = 1:N
        plot(px(i,k), py(i,k), 'o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', 'k', ...
           'MarkerSize', mk_size, 'Color', [cols(i,:), alpha_val]);
    end
end
end

%% ======================================================
%  Formation animation (Figure 2)
%% ======================================================
function animate_formation(res1, res2, t, N, cols, soc_min, L_topo)

T_sim    = length(t);
skip     = 4;            % animate every 4 steps
trail    = 120;          % trail length [steps]
win      = 6.0;          % camera half-width [m]

% Find communication graph edges (i < j, undirected)
net_edges = [];
for i = 1:N
    for j = i+1:N
        if L_topo(i,j) ~= 0
            net_edges = [net_edges; i, j]; %#ok<AGROW>
        end
    end
end
n_net = size(net_edges, 1);

fig = figure(2);
fig.Name     = 'Formation Animation';
fig.Position = [80 80 1400 640];
clf;

ax      = gobjects(1, 2);
res_all = {res1, res2};
labels  = {'Proposed (Energy-Aware)', 'Baseline'};

% ---- Pre-create graphics handles ----
h_net    = cell(2, n_net);  % communication network edges (gray dashed)
h_trail  = cell(2, N);
h_edge   = cell(2, N);      % geometric pentagon outline (black solid)
h_socbar = cell(2, N);
h_agent  = cell(2, N);
h_soclbl = cell(2, N);

for p = 1:2
    ax(p) = subplot(1, 2, p);
    hold(ax(p), 'on');
    axis(ax(p), 'equal');
    grid(ax(p), 'on');
    xlabel(ax(p), 'x [m]');
    ylabel(ax(p), 'y [m]');

    % 1) Network edges — light gray dashed (drawn first = behind everything)
    for e = 1:n_net
        h_net{p,e} = plot(ax(p), nan(1,2), nan(1,2), '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 1.0);
    end

    % 2) Agent trails
    for i = 1:N
        h_trail{p,i} = plot(ax(p), nan, nan, '-', 'Color', [cols(i,:), 0.25], 'LineWidth', 1.2);
    end

    % 3) Geometric pentagon outline (1→2→3→4→5→1, black solid)
    for e = 1:N
        h_edge{p,e} = plot(ax(p), nan(1,2), nan(1,2), 'k-', 'LineWidth', 1.8);
    end

    % 4) SOC bar + agent dot + SOC label
    for i = 1:N
        h_socbar{p,i} = fill(ax(p), nan(1,4), nan(1,4), ...
            soc_color(res_all{p}.soc(i,1)), ...
            'EdgeColor', 'k', 'LineWidth', 0.8, 'FaceAlpha', 0.85);
        h_agent{p,i} = plot(ax(p), nan, nan, 'o', ...
            'MarkerSize', 14, ...
            'MarkerFaceColor', cols(i,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        h_soclbl{p,i} = text(ax(p), nan, nan, '', ...
            'FontSize', 7.5, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

% ---- Animation loop ----
for k = 1:skip:T_sim
    for p = 1:2
        res = res_all{p};

        % Moving camera centred on formation
        cx = mean(res.px(:,k));
        cy = mean(res.py(:,k));
        xlim(ax(p), [cx-win, cx+win]);
        ylim(ax(p), [cy-win, cy+win]);

        k0 = max(1, k - trail);

        % -- Network edges (gray dashed) --
        for e = 1:n_net
            ii = net_edges(e,1);
            jj = net_edges(e,2);
            set(h_net{p,e}, ...
                'XData', [res.px(ii,k), res.px(jj,k)], ...
                'YData', [res.py(ii,k), res.py(jj,k)]);
        end

        % -- Per-agent: trail, SOC bar, dot, label --
        for i = 1:N
            px_i = res.px(i,:);
            py_i = res.py(i,:);
            s    = res.soc(i,k);

            % Trail
            set(h_trail{p,i}, 'XData', px_i(k0:k), 'YData', py_i(k0:k));

            % SOC bar (filled rectangle below agent, width ∝ SOC)
            bw = 0.5; bh = 0.14;
            bx = px_i(k) - bw/2;
            by = py_i(k) - 1.1;
            set(h_socbar{p,i}, ...
                'XData', [bx, bx+bw*s, bx+bw*s, bx, bx], ...
                'YData', [by, by, by+bh, by+bh, by], ...
                'FaceColor', soc_color(s));

            % Agent dot
            set(h_agent{p,i}, 'XData', px_i(k), 'YData', py_i(k));

            % SOC label
            set(h_soclbl{p,i}, ...
                'Position', [px_i(k), py_i(k)-0.72, 0], ...
                'String', sprintf('%d%%', round(s*100)));
        end

        % -- Geometric pentagon outline (1→2→3→4→5→1) --
        for e = 1:N
            j = mod(e, N) + 1;
            set(h_edge{p,e}, ...
                'XData', [res.px(e,k), res.px(j,k)], ...
                'YData', [res.py(e,k), res.py(j,k)]);
        end

        title(ax(p), sprintf('%s   t = %.2f s', labels{p}, t(k)), ...
            'FontSize', 11, 'FontWeight', 'bold');
    end

    sgtitle(fig, 'Pentagon Formation — Energy-Aware vs. Baseline', ...
        'FontSize', 13, 'FontWeight', 'bold');

    drawnow limitrate;
end

end

%% ======================================================
%  SOC → color helper  (green=high, red=low)
%% ======================================================
function c = soc_color(s)
    % s in [0,1] → RGB: green (1) to red (0)
    s = max(0, min(1, s));
    c = [1-s, s, 0.1];
end

%% ======================================================
%  Internal simulation runner
%% ======================================================
function res = run_sim(A, B, C, E, F, L_topo, K, L_obs, delta, x0, xhat0, soc0, W_noise, V_noise, T_sim, alpha, eta_ctrl, eta_prop, soc_min, soc_max, method)

    N    = size(L_topo, 1);
    n    = size(A, 1);
    m    = size(B, 2);
    q    = size(C, 1);
    
    x    = x0;
    xhat = xhat0;
    soc  = soc0;
    
    % Fixed weight for Baseline
    w_upper = exp(-2 * alpha / soc_max);
    
    % Pre-allocate
    form_err   = zeros(1, T_sim);
    u_total    = zeros(1, T_sim);
    cum_energy = zeros(1, T_sim);
    soc_log    = zeros(N, T_sim);
    px_log     = zeros(N, T_sim);
    py_log     = zeros(N, T_sim);
    u_each     = zeros(N, T_sim);   % per-agent ||u_i(k)||
    v_each     = zeros(N, T_sim);   % per-agent speed ||v_i(k)||
    cum_e      = 0;
    
    for k = 1:T_sim
    
        %-- Shifted (formation error) states --
        % z_i = x_i - delta_i: if all z_i coincide → pentagon achieved
        zhat = xhat - delta;    % n x N
    
        %-- Build edge weight matrix (N x N, positive) --
        W_mat = zeros(N, N);
        if strcmp(method, 'proposed')
            soc_k = max(soc, soc_min);
            for i = 1:N
                for j = 1:N
                    if i ~= j && L_topo(i,j) ~= 0
                        W_mat(i,j) = exp(-alpha * (1/soc_k(i) + 1/soc_k(j)));
                    end
                end
            end
        else
            for i = 1:N
                for j = 1:N
                    if i ~= j && L_topo(i,j) ~= 0
                        W_mat(i,j) = w_upper;
                    end
                end
            end
        end
    
        %-- Control: u_i = K * sum_j w_ij*(zhat_i - zhat_j) --
        u = zeros(m, N);
        for i = 1:N
            weighted_sum = zeros(n, 1);
            for j = 1:N
                if W_mat(i,j) > 0
                    weighted_sum = weighted_sum + W_mat(i,j) * (zhat(:,i) - zhat(:,j));
                end
            end
            u(:,i) = K * weighted_sum;
        end
    
        %-- Formation error: deviation of z_i from mean z --
        z_mean        = mean(x - delta, 2);
        form_err(k)   = sum(vecnorm(x - delta - repmat(z_mean, 1, N), 2, 1));
    
        %-- Log positions --
        px_log(:,k)   = x(1,:)';
        py_log(:,k)   = x(3,:)';
    
        %-- Control effort metrics --
        u_norms       = vecnorm(u, 2, 1);
        u_total(k)    = sum(u_norms);
        cum_e         = cum_e + sum(u_norms.^2);
        cum_energy(k) = cum_e;
        u_each(:,k)   = u_norms(:);
    
        %-- SOC update --
        %   Control cost:   eta_ctrl * ||u_i||^2
        %   Propulsion cost: eta_prop * (vx_i^2 + vy_i^2)  [agents keep moving]
        v_sq          = x(2,:).^2 + x(4,:).^2;    % ||velocity||^2 for each agent
        v_each(:,k)   = sqrt(v_sq(:));
        soc           = soc - eta_ctrl * u_norms(:).^2 - eta_prop * v_sq(:);
        soc           = max(soc, 0);
        soc_log(:,k)  = soc;
    
        %-- Step 2: Measure y(k) from CURRENT x(k) [before plant update] --
        y_k           = C*x + F*squeeze(V_noise(:,:,k));
    
        %-- Step 3: Observer update --
        for i = 1:N
            xhat(:,i) = A*xhat(:,i) + B*u(:,i) + L_obs * (y_k(:,i) - C*xhat(:,i));
        end
    
        %-- Step 4: Plant update --
        x = A*x + B*u + E*squeeze(W_noise(:, :, k));
    
    end
    
    res.form_err   = form_err;
    res.u_total    = u_total;
    res.cum_energy = cum_energy;
    res.soc        = soc_log;
    res.px         = px_log;
    res.py         = py_log;
    res.u_each     = u_each;
    res.v_each     = v_each;

end
