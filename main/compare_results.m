function compare_results(varargin)
% COMPARE_RESULTS  Generate all manuscript figures from saved simulation results.
%
%   compare_results()             -- generates figures for all three scenarios
%   compare_results('scenario','B_rapid')  -- single scenario
%
%   Requires run_dyas_sim() to have been run first.
%   Figures saved to results/ as .png and .pdf

p = inputParser;
addParameter(p, 'scenario', 'all');
parse(p, varargin{:});
sc_arg = p.Results.scenario;

rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(rootdir));

if strcmp(sc_arg, 'all')
    sc_list = {'A_mild', 'B_rapid', 'C_prefatigued'};
else
    sc_list = {sc_arg};
end

% Load all results
all_res = struct();
for s = 1:numel(sc_list)
    sc = sc_list{s};
    f  = fullfile(rootdir, 'results', ['results_', sc, '.mat']);
    if ~exist(f, 'file')
        error('Results file not found: %s\nRun run_dyas_sim() first.', f);
    end
    tmp = load(f);
    all_res.(sc).dyas = tmp.results_dyas;
    all_res.(sc).full = tmp.results_full;
    all_res.(sc).S    = tmp.S;
    all_res.(sc).P    = tmp.P;
end

P = all_res.(sc_list{1}).P;
colors = struct('full', [0.3 0.3 0.8], 'dyas', [0.8 0.2 0.2]);
sc_labels = {'A – Mild', 'B – Rapid', 'C – Pre-fat.'};

% =========================================================================
%  Figure 1: Per-cycle mean FES capacity  (FigL1)
% =========================================================================
figure('Name','FigL1_fatigue_cycles','Position',[100 100 900 300]);
tiledlayout(1, numel(sc_list), 'TileSpacing','compact','Padding','compact');

for s = 1:numel(sc_list)
    sc   = sc_list{s};
    r_d  = all_res.(sc).dyas;
    r_f  = all_res.(sc).full;

    % Mean over FES channels only (channels 4:8)
    phi_mean_d = mean(r_d.phi_per_cycle(P.nm+1:end, :), 1);
    phi_mean_f = mean(r_f.phi_per_cycle(P.nm+1:end, :), 1);

    nexttile;
    cycles = 1:P.N_cycles;
    plot(cycles, phi_mean_f, '--', 'Color', colors.full, 'LineWidth', 1.5);
    hold on;
    plot(cycles, phi_mean_d, '-',  'Color', colors.dyas, 'LineWidth', 1.5);
    xlabel('Cycle'); ylabel('\bar{\phi}');
    title(sc_labels{s});
    ylim([0 1]); grid on;
    if s == 1
        legend('Full (n_w=8)', 'DyAS (l^*=3)', 'Location','southwest');
    end
    % Annotate gain
    gain_pp = (phi_mean_d(end) - phi_mean_f(end)) * 100;
    text(P.N_cycles*0.6, 0.1, sprintf('+%.1f pp', gain_pp), ...
         'Color', colors.dyas, 'FontSize', 9, 'FontWeight','bold');
end
sgtitle('Per-cycle FES endurance: Full vs DyAS allocator (20-cycle STS)');
save_figure(gcf, rootdir, 'FigL1_fatigue_cycles');

% =========================================================================
%  Figure 2: Channel-wise fatigue timeseries for Case B  (FigL2)
% =========================================================================
sc  = 'B_rapid';
r_d = all_res.(sc).dyas;
r_f = all_res.(sc).full;

fes_labels = {'TA','SOL','QUAD','HAMS','GLUT'};
cmap = lines(P.nf);

figure('Name','FigL2_fatigue_timeseries','Position',[100 100 1000 350]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for panel = 1:2
    nexttile;
    if panel == 1
        phi_data = r_f.phi(P.nm+1:end, :);
        ttl = 'Full (n_w=8)';
    else
        phi_data = r_d.phi(P.nm+1:end, :);
        ttl = 'DyAS (l^*=3)';
    end
    t_vec = r_d.t;
    for ch = 1:P.nf
        plot(t_vec, phi_data(ch,:), 'Color', cmap(ch,:), 'LineWidth', 1.3);
        hold on;
    end
    xlabel('Time (s)'); ylabel('\phi_i');
    title(ttl);
    ylim([0.1 1.05]); grid on;
    legend(fes_labels, 'Location','southwest', 'FontSize',7);
end
sgtitle('Case B (Rapid) – Individual FES channel fatigue over 140 s');
save_figure(gcf, rootdir, 'FigL2_fatigue_timeseries');

% =========================================================================
%  Figure 3: Per-cycle joint RMSE  (FigL3)
% =========================================================================
joint_labels = {'Ankle','Knee','Hip'};
line_style_d = {'-','-.',':'};
line_style_f = {'--','-..',':'};

figure('Name','FigL3_rmse_per_cycle','Position',[100 100 900 300]);
tiledlayout(1, numel(sc_list), 'TileSpacing','compact','Padding','compact');

for s = 1:numel(sc_list)
    sc  = sc_list{s};
    r_d = all_res.(sc).dyas;
    r_f = all_res.(sc).full;
    cycles = 1:P.N_cycles;

    nexttile; hold on;
    for j = 1:P.n
        plot(cycles, r_f.rmse_per_cycle(j,:)*180/pi, '--', ...
             'Color', colors.full * (0.5 + 0.1*j), 'LineWidth', 1.2);
        plot(cycles, r_d.rmse_per_cycle(j,:)*180/pi, '-', ...
             'Color', colors.dyas * (0.5 + 0.1*j), 'LineWidth', 1.2);
    end
    xlabel('Cycle'); ylabel('RMSE (deg)');
    title(sc_labels{s}); grid on;
    if s == 1
        legend({'Ankle Full','Ankle DyAS','Knee Full','Knee DyAS',...
                'Hip Full','Hip DyAS'}, 'FontSize',6, 'Location','northeast');
    end
end
sgtitle('Per-cycle joint tracking RMSE over 20 STS cycles');
save_figure(gcf, rootdir, 'FigL3_rmse_per_cycle');

% =========================================================================
%  Figure 4: Eigenvalue evolution  (FigL4) — Case B
% =========================================================================
sc  = 'B_rapid';
r_d = all_res.(sc).dyas;

% Re-simulate to capture eigenvalue history
[eig_history, l_history] = get_eigenvalue_history(r_d, P, ...
    dyas_params(), fatigue_scenarios(sc));

figure('Name','FigL4_eigenvalue_evolution','Position',[100 100 700 450]);

subplot(2,1,1);
t_dyas = (0:size(eig_history,2)-1) * dyas_params().dt_dyas;
norm_eig = eig_history ./ sum(eig_history, 1);
for j = 1:min(P.nw, 6)
    if j <= 3
        lw = 2.0; ls = '-';
    else
        lw = 1.0; ls = '--';
    end
    plot(t_dyas, norm_eig(j,:), ls, 'LineWidth', lw);
    hold on;
end
ylabel('\lambda_j / tr(C)');
title('Case B (Rapid) – GCM eigenvalue evolution');
legend({'\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','\lambda_6'}, ...
       'Location','northeast','FontSize',8);
grid on; ylim([0 1]);
yline(0.05, 'k:', 'LineWidth', 1);

subplot(2,1,2);
t_full = r_d.t;
stairs(t_full, r_d.l_star, 'r-', 'LineWidth', 2);
hold on;
yline(P.nw, 'Color',[0.5 0.5 0.5], 'LineStyle','--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('l^*');
yticks([1 2 3 4 5 6 7 8]);
ylim([0.5 P.nw+0.5]);
legend('DyAS l^*', sprintf('Full n_w=%d', P.nw), 'Location','east');
title('Active subspace dimension vs time');
grid on;

save_figure(gcf, rootdir, 'FigL4_eigenvalue_evolution');

% =========================================================================
%  Figure 5: End-of-session FES capacity bar chart  (FigL5)
% =========================================================================
phi_end_full = zeros(1, numel(sc_list));
phi_end_dyas = zeros(1, numel(sc_list));

for s = 1:numel(sc_list)
    sc = sc_list{s};
    r_d = all_res.(sc).dyas;
    r_f = all_res.(sc).full;
    phi_end_full(s) = mean(r_f.phi_per_cycle(P.nm+1:end, end));
    phi_end_dyas(s) = mean(r_d.phi_per_cycle(P.nm+1:end, end));
end

figure('Name','FigL5_endurance_bar','Position',[100 100 500 350]);
x  = 1:numel(sc_list);
bw = 0.35;
b1 = bar(x - bw/2, phi_end_full, bw, 'FaceColor', colors.full);
hold on;
b2 = bar(x + bw/2, phi_end_dyas, bw, 'FaceColor', colors.dyas);

% Annotate gains
for s = 1:numel(sc_list)
    gain = (phi_end_dyas(s) - phi_end_full(s)) * 100;
    text(s + bw/2, phi_end_dyas(s) + 0.02, sprintf('+%.1f pp', gain), ...
         'HorizontalAlignment','center', 'FontSize', 9, ...
         'Color', colors.dyas, 'FontWeight','bold');
end

set(gca, 'XTick', x, 'XTickLabel', sc_labels);
ylabel('Mean FES capacity \bar{\phi} at cycle 20');
legend([b1 b2], {'Full (n_w=8)', 'DyAS (l^*=3)'}, 'Location','north');
title('FES endurance at session end (20-cycle STS)');
ylim([0 1.1]); grid on;
save_figure(gcf, rootdir, 'FigL5_endurance_bar');

fprintf('\nAll figures saved to results/\n');
end


% =========================================================================
%  Helper: extract eigenvalue history from stored phi trajectory
% =========================================================================
function [eig_hist, l_hist] = get_eigenvalue_history(R, P, D, S)
% Re-run the DyAS update on stored phi and tau trajectories to extract
% eigenvalue evolution (since the main loop doesn't store GCM snapshots).

n_steps = size(R.phi, 2);
n_upd   = floor(P.T_session / D.dt_dyas);
eig_hist = zeros(P.nw, n_upd);
l_hist   = zeros(1, n_upd);

C_gcm = eye(P.nw);
Pi    = eye(P.nw);
upd   = 0;

for k = 1:n_steps
    t = (k-1)*P.dt;
    if t >= upd * D.dt_dyas && upd < n_upd
        upd = upd + 1;
        phi_k = R.phi(:,k);
        Phi_k = diag(phi_k);
        tau_k = R.tau(:,k);

        [~, l_star, C_gcm, Pi] = dyas_update(Phi_k, tau_k, C_gcm, Pi, P, D);

        [~, lam] = eig(C_gcm, 'vector');
        lam = sort(lam, 'descend');
        lam = max(lam, 0);
        eig_hist(:, upd) = lam;
        l_hist(upd) = l_star;
    end
end
end


% =========================================================================
%  Helper: save figure as PNG and PDF
% =========================================================================
function save_figure(fig, rootdir, fname)
outdir = fullfile(rootdir, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
print(fig, fullfile(outdir, fname), '-dpng', '-r150');
print(fig, fullfile(outdir, fname), '-dpdf');
fprintf('  Saved figure: %s\n', fname);
end
