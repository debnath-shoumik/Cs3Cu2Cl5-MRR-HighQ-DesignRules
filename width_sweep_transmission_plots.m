% %% Load and process your data
% clear; clc; close all;
% 
% % Load the data
% load('ring_width_sweep.mat');
% 
% base_width_m = 500e-09;  % 500 nm in meters
% actual_widths_nm = (base_width_m + width_vals) * 1e9;  % Convert to nm
% 
% % Convert wavelength to nm for easier reading
% lambda_nm = lambda * 1e9;  % meters to nanometers
% 
% 
% 
% %% Individual plots for each width
% figure('Position', [50, 50, 1400, 800]);
% for i = 1:length(actual_widths_nm)
%     subplot(2, 3, i);
% 
%     % Plot both drop and through
%     plot(lambda_nm, T_drop(i,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Drop Port');
%     hold on;
%     plot(lambda_nm, T_through(i,:), 'r-', 'LineWidth', 2, 'DisplayName', 'Through Port');
%     hold off;
% 
%     xlabel('Wavelength (nm)', 'FontSize', 10, 'FontWeight', 'bold');
%     ylabel('Transmission', 'FontSize', 10, 'FontWeight', 'bold');
%     title(sprintf('Width = %.0f nm', actual_widths_nm(i)), 'FontSize', 12, 'FontWeight', 'bold');
%     legend('Location', 'best', 'FontSize', 9);
%     grid on;
%     xlim([min(lambda_nm), max(lambda_nm)]);
%     ylim([0, 1]);
% end
% 
% 
% %% Representative spectrum (Width = 600 nm)
% 
% idx = 2;  % 600 nm case
% 
% figure('Color','w','Position',[300 300 700 500]);
% 
% plot(lambda_nm, T_drop(idx,:), ...
%     'Color',[0 0.2 0.6], ...
%     'LineWidth',1.8);
% hold on;
% 
% plot(lambda_nm, T_through(idx,:), ...
%     'Color',[0.6 0 0], ...
%     'LineWidth',1.8);
% 
% xlabel('Wavelength (nm)', ...
%     'FontName','Times New Roman','FontSize',14);
% 
% ylabel('Transmission', ...
%     'FontName','Times New Roman','FontSize',14);
% 
% set(gca,'FontName','Times New Roman', ...
%     'FontSize',13, ...
%     'LineWidth',1.2);
% 
% box off;
% legend({'Drop Port','Through Port'}, ...
%     'Location','northeast', ...
%     'FontSize',12);
% 
% xlim([min(lambda_nm) max(lambda_nm)]);
% ylim([0 1]);
% 
% 
% %% Width comparison – Drop Port Only
% 
% figure('Color','w','Position',[300 300 700 500]);
% hold on;
% 
% colors = lines(length(actual_widths_nm));
% 
% for i = 1:length(actual_widths_nm)
%     plot(lambda_nm, T_drop(i,:), ...
%         'LineWidth',1.6, ...
%         'Color',colors(i,:));
% end
% 
% xlabel('Wavelength (nm)', ...
%     'FontName','Times New Roman','FontSize',14);
% 
% ylabel('Drop Port Transmission', ...
%     'FontName','Times New Roman','FontSize',14);
% 
% set(gca,'FontName','Times New Roman', ...
%     'FontSize',13, ...
%     'LineWidth',1.2);
% 
% box off;
% 
% legend(arrayfun(@(x) sprintf('%d nm',round(x)), ...
%     actual_widths_nm, ...
%     'UniformOutput',false), ...
%     'Location','northeast', ...
%     'FontSize',11);
% 
% xlim([min(lambda_nm) max(lambda_nm)]);
% ylim([0 1]);




%% ================================================================
%  FIGURE 2 — Width Sweep Transmission Spectra
%  Cs₃Cu₂Cl₅ Microring Resonator | Optical Materials Express
%  Two-panel: (a) representative spectra, (b) drop port comparison
%  ================================================================
clear; clc; close all;

load('ring_width_sweep.mat');

base_width_m     = 500e-09;
actual_widths_nm = (base_width_m + width_vals) * 1e9;
lambda_nm        = lambda * 1e9;

%% ——— Color palette ———
C = lines(7);
%  1: blue  2: orange  3: yellow  4: purple  5: green  6: cyan  7: maroon

% Panel (a): 3 representative widths — through = lighter/thinner, drop = bold
rep_idx    = [1, 2, 4];           % 500, 600, 800 nm
rep_labels = {'500', '600', '800'};
rep_colors = [C(1,:); C(2,:); C(7,:)];

% Panel (b): all 5 widths
all_colors = [C(1,:); C(2,:); C(4,:); C(5,:); C(7,:)];

%% ——— Figure ———
fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 16 17]);

% =========================
%  Panel (a)
% =========================
ax1 = subplot(2,1,1);
hold on;

h = gobjects(6,1);
for k = 1:3
    i = rep_idx(k);
    % Through: same color but semi-transparent (alpha via 4th color channel trick)
    % Since MATLAB line alpha is limited, use a lighter tint instead
    tint = 1 - 0.45*(1 - rep_colors(k,:));   % 45 % toward white
    h(2*k-1) = plot(lambda_nm, T_through(i,:), '-', ...
        'Color', tint, 'LineWidth', 1.1);
    h(2*k)   = plot(lambda_nm, T_drop(i,:), '-', ...
        'Color', rep_colors(k,:), 'LineWidth', 1.8);
end
hold off;

xlim([min(lambda_nm) max(lambda_nm)]);
ylim([0 1.05]);
ylabel('Transmission');
set(ax1, 'XTickLabel', []);

% Legend
leg_entries = {};
h_leg = [];
for k = 1:3
    leg_entries{end+1} = sprintf('%s nm (through)', rep_labels{k});
    leg_entries{end+1} = sprintf('%s nm (drop)',    rep_labels{k});
    h_leg = [h_leg; h(2*k-1); h(2*k)];
end
legend(h_leg, leg_entries, ...
    'Location','east', 'FontSize',12, 'NumColumns',1, 'Box','off');

text(0.025, 0.96, '(a)', 'Units','normalized', ...
    'FontSize',21, 'FontWeight','bold', ...
    'FontName','Times New Roman', 'VerticalAlignment','top');

set(ax1, 'FontName','Times New Roman', 'FontSize',21, ...
    'TickDir','in', 'LineWidth',0.9, ...
    'XMinorTick','on', 'YMinorTick','on', 'Box','on');

% =========================
%  Panel (b)
% =========================
ax2 = subplot(2,1,2);
hold on;
for i = 1:length(actual_widths_nm)
    plot(lambda_nm, T_drop(i,:), '-', ...
        'Color', all_colors(i,:), 'LineWidth', 1.5);
end
hold off;

xlim([min(lambda_nm) max(lambda_nm)]);
ylim([0 1.05]);
xlabel('Wavelength (nm)');
ylabel('Drop port transmission');

leg_str = arrayfun(@(x) sprintf('W_r = %d nm', round(x)), ...
    actual_widths_nm, 'UniformOutput', false);
legend(leg_str, 'Location','northeast', 'FontSize',12, 'Box','off');

text(0.025, 0.96, '(b)', 'Units','normalized', ...
    'FontSize',21, 'FontWeight','bold', ...
    'FontName','Times New Roman', 'VerticalAlignment','top');

set(ax2, 'FontName','Times New Roman', 'FontSize',21, ...
    'TickDir','in', 'LineWidth',0.9, ...
    'XMinorTick','on', 'YMinorTick','on', 'Box','on');

% ——— Layout ———
ax1.Position = [0.13 0.55 0.82 0.41];
ax2.Position = [0.13 0.08 0.82 0.41];

% ——— Export ———
exportgraphics(fig, 'Fig2_width_sweep.png',  'Resolution', 600);
exportgraphics(fig, 'Fig2_width_sweep.pdf',  'ContentType', 'vector');
fprintf('Done — Fig2_width_sweep.png / .pdf\n');