% %% Load and process your data
% clear; clc; close all;
% 
% % Load the data
% load('ring_width_sweep.mat');
% 
% % Check what width_vals contains
% fprintf('width_vals contains (meters):\n');
% disp(width_vals);
% 
% % Convert to actual widths in nm
% base_width_m = 500e-09;  % 500 nm in meters
% actual_widths_nm = (base_width_m + width_vals) * 1e9;  % Convert to nm
% 
% fprintf('\nActual widths in nm:\n');
% disp(actual_widths_nm);
% 
% % Convert wavelength to nm for easier reading
% lambda_nm = lambda * 1e9;  % meters to nanometers
% 
% %% Display parameters FOR EACH WIDTH
% fprintf('\n==============================================\n');
% fprintf('MICRO RING RESONATOR ANALYSIS\n');
% fprintf('==============================================\n');
% 
% % Pre-allocate results structure
% all_results = struct();
% 
% % Loop through each width and calculate parameters
% for i = 1:length(actual_widths_nm)
%     fprintf('\n--- Width = %.0f nm ---\n', actual_widths_nm(i));
% 
%     % Get spectra for this width
%     drop_spectrum = T_drop(i,:);
%     through_spectrum = T_through(i,:);
% 
%     % 1. Find resonance wavelength (min in drop port)
%     [min_drop, idx_min] = min(drop_spectrum);
%     resonance_lambda = lambda_nm(idx_min);
% 
%     % 2. Find FSR (find adjacent resonances in through port)
%     % Use findpeaks to locate resonance peaks in through port
%     [pks, locs] = findpeaks(through_spectrum, 'MinPeakProminence', 0.1);
% 
%     if length(locs) >= 2
%         FSR_values = diff(lambda_nm(locs));
%         avg_FSR = mean(abs(FSR_values));
%         fprintf('Found %d resonances, FSR = %.3f nm\n', length(locs), avg_FSR);
%     else
%         avg_FSR = NaN;
%         fprintf('Only found %d resonance(s)\n', length(locs));
%     end
% 
%     % 3. Calculate FWHM PROPERLY with interpolation
%     % Find the half-max value
%     max_val = max(drop_spectrum);
%     min_val = min_drop;
%     half_max = (max_val + min_val) / 2;
% 
%     % Find where the spectrum crosses half-max
%     % Interpolate for accurate crossing points
%     fwhm_nm = calculate_fwhm(lambda_nm, drop_spectrum, half_max);
% 
%     % 4. Calculate Q-factor
%     if ~isnan(fwhm_nm) && fwhm_nm > 0
%         Q_factor = resonance_lambda / fwhm_nm;
%     else
%         Q_factor = NaN;
%     end
% 
%     % 5. Drop efficiency (1 - min transmission at resonance)
%     drop_efficiency = 1 - min_drop;
% 
%     % 6. Insertion Loss (at through port maximum, away from resonance)
%     % Find maximum away from resonance (e.g., at edges)
%     edge_samples = 10;
%     max_through = max(through_spectrum([1:edge_samples, end-edge_samples:end]));
%     insertion_loss_db = -10 * log10(max_through);
% 
%     % 7. Extinction Ratios
%     ER_through_db = -10 * log10(min(through_spectrum)/max(through_spectrum));
%     ER_drop_db = -10 * log10(min_val/max_val);
% 
%     % 8. Finesse
%     if ~isnan(avg_FSR) && ~isnan(fwhm_nm) && fwhm_nm > 0
%         finesse = avg_FSR / fwhm_nm;
%     else
%         finesse = NaN;
%     end
% 
%     % 9. Peak Contrast
%     peak_contrast = (max_through - min(through_spectrum)) / max_through;
% 
%     % Display results
%     fprintf('Resonance wavelength  = %.3f nm\n', resonance_lambda);
%     fprintf('Average FSR           = %.3f nm\n', avg_FSR);
%     fprintf('FWHM                  = %.4f nm\n', fwhm_nm);
%     fprintf('Q-factor              = %.1f\n', Q_factor);
%     fprintf('Drop contrast         = %.3f\n', drop_efficiency);
%     fprintf('Insertion Loss        = %.2f dB\n', insertion_loss_db);
%     fprintf('ER (Through)          = %.2f dB\n', ER_through_db);
%     fprintf('ER (Drop)             = %.2f dB\n', ER_drop_db);
%     fprintf('Finesse               = %.2f\n', finesse);
%     fprintf('Peak Contrast         = %.3f\n\n', peak_contrast);
% 
%     % Store results
%     all_results(i).Width_nm = actual_widths_nm(i);
%     all_results(i).Resonance_nm = resonance_lambda;
%     all_results(i).FSR_nm = avg_FSR;
%     all_results(i).FWHM_nm = fwhm_nm;
%     all_results(i).Q_factor = Q_factor;
%     all_results(i).DropContrast = drop_efficiency;
%     all_results(i).InsertionLoss_dB = insertion_loss_db;
%     all_results(i).ER_Through_dB = ER_through_db;
%     all_results(i).ER_Drop_dB = ER_drop_db;
%     all_results(i).Finesse = finesse;
%     all_results(i).PeakContrast = peak_contrast;
%     all_results(i).DropSpectrum = drop_spectrum;
%     all_results(i).ThroughSpectrum = through_spectrum;
% end
% 
% %% Create and display results table
% fprintf('\n==============================================\n');
% fprintf('SUMMARY TABLE\n');
% fprintf('==============================================\n');
% 
% % Convert structure to table for nice display
% results_table = struct2table(all_results);
% disp(results_table(:,1:10)); % Display first 10 columns
% 
% %% Plot trends vs width
% % figure('Position', [100, 100, 1200, 800]);
% % 
% % % Subplot 1: Resonance wavelength vs width
% % subplot(2,3,1);
% % resonances = [all_results.Resonance_nm];
% % plot(actual_widths_nm, resonances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('Resonance λ (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % title('Resonance Wavelength vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% % 
% % % Subplot 2: Q-factor vs width
% % subplot(2,3,2);
% % Q_factors = [all_results.Q_factor];
% % plot(actual_widths_nm, Q_factors, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('Q-factor', 'FontSize', 11, 'FontWeight', 'bold');
% % title('Q-factor vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% % 
% % % Subplot 3: FSR vs width
% % subplot(2,3,3);
% % FSRs = [all_results.FSR_nm];
% % plot(actual_widths_nm, FSRs, 'go-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('FSR (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % title('FSR vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% % 
% % % Subplot 4: ER Drop vs width
% % subplot(2,3,4);
% % ER_drops = [all_results.ER_Drop_dB];
% % plot(actual_widths_nm, ER_drops, 'mo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('ER Drop (dB)', 'FontSize', 11, 'FontWeight', 'bold');
% % title('Drop Port ER vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% % 
% % % Subplot 5: Finesse vs width
% % subplot(2,3,5);
% % Finesses = [all_results.Finesse];
% % plot(actual_widths_nm, Finesses, 'co-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'c');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('Finesse', 'FontSize', 11, 'FontWeight', 'bold');
% % title('Finesse vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% % 
% % % Subplot 6: FWHM vs width
% % subplot(2,3,6);
% % FWHMs = [all_results.FWHM_nm];
% % plot(actual_widths_nm, FWHMs, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
% % xlabel('Ring Width (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % ylabel('FWHM (nm)', 'FontSize', 11, 'FontWeight', 'bold');
% % title('FWHM vs Width', 'FontSize', 12, 'FontWeight', 'bold');
% % grid on;
% 
% 
% 
% %% Helper function for accurate FWHM calculation
% function fwhm = calculate_fwhm(x, y, half_max)
%     % Find all indices where y crosses half_max
%     cross_idx = find(diff(sign(y - half_max)) ~= 0);
% 
%     if length(cross_idx) >= 2
%         % Interpolate for accurate crossing points
%         x1 = interp1(y(cross_idx(1):cross_idx(1)+1), ...
%                      x(cross_idx(1):cross_idx(1)+1), half_max);
%         x2 = interp1(y(cross_idx(2):cross_idx(2)+1), ...
%                      x(cross_idx(2):cross_idx(2)+1), half_max);
%         fwhm = abs(x2 - x1);
%     else
%         % Try alternative method: find closest points to half_max
%         [~, idx1] = min(abs(y(1:end-1) - half_max));
%         [~, idx2] = min(abs(y(2:end) - half_max));
%         fwhm = abs(x(idx2+1) - x(idx1));
%     end
% 
%     % Sanity check
%     if fwhm > 100 || fwhm < 0.001
%         fwhm = NaN;
%     end
% end
% 
% %% Save results to Excel file
% output_filename = 'ring_width_sweep_results.xlsx';
% writetable(results_table, output_filename);
% fprintf('\nResults saved to: %s\n', output_filename);
% 
% 
% 
% %% ===================== SEPARATE FULL FIGURES =====================
% 
% skyblue = [0 0.4470 0.7410];   % MATLAB default sky-blue
% 
% % 1. Resonance wavelength vs width
% y = [all_results.Resonance_nm];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Resonance wavelength (nm)');
% title('Resonance wavelength vs ring width');
% ylim([0.98*min(y) 1.02*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 2. Q-factor vs width
% y = [all_results.Q_factor];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Q factor');
% title('Q factor vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 3. FSR vs width
% y = [all_results.FSR_nm];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('FSR (nm)');
% title('FSR vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 4. Drop extinction ratio vs width
% y = [all_results.ER_Drop_dB];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Drop extinction ratio (dB)');
% title('Drop port extinction ratio vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 5. Finesse vs width
% y = [all_results.Finesse];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Finesse');
% title('Finesse vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 6. FWHM vs width
% y = [all_results.FWHM_nm];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('FWHM (nm)');
% title('FWHM vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 7. Insertion loss vs width
% y = [all_results.InsertionLoss_dB];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Insertion loss (dB)');
% title('Insertion loss vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;
% 
% % 8. Peak contrast vs width
% y = [all_results.PeakContrast];
% figure;
% plot(actual_widths_nm, y, 'o-', ...
%     'Color', skyblue, 'MarkerFaceColor', skyblue);
% xlabel('Ring width (nm)');
% ylabel('Peak contrast');
% title('Peak contrast vs ring width');
% ylim([0.9*min(y) 1.1*max(y)]);
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid off;






%% ================================================================
%  FIGURE 3 — Extracted Parameters vs Ring Width
%  Cs₃Cu₂Cl₅ Microring Resonator | Optical Materials Express
%  Three-panel: (a) Q-factor, (b) Drop ER, (c) Finesse
%  ================================================================
clear; clc; close all;

load('ring_width_sweep.mat');

base_width_m     = 500e-09;
actual_widths_nm = (base_width_m + width_vals) * 1e9;
lambda_nm        = lambda * 1e9;

%% ——— Compute parameters (from your analysis script) ———
all_results = struct();
for i = 1:length(actual_widths_nm)
    drop_spectrum    = T_drop(i,:);
    through_spectrum = T_through(i,:);

    [min_drop, idx_min] = min(drop_spectrum);
    resonance_lambda = lambda_nm(idx_min);

    [~, locs] = findpeaks(through_spectrum, 'MinPeakProminence', 0.1);
    if length(locs) >= 2
        avg_FSR = mean(abs(diff(lambda_nm(locs))));
    else
        avg_FSR = NaN;
    end

    max_val  = max(drop_spectrum);
    half_max = (max_val + min_drop) / 2;
    fwhm_nm  = calculate_fwhm(lambda_nm, drop_spectrum, half_max);

    if ~isnan(fwhm_nm) && fwhm_nm > 0
        Q_factor = resonance_lambda / fwhm_nm;
    else
        Q_factor = NaN;
    end

    ER_drop_db = -10 * log10(min_drop / max_val);

    if ~isnan(avg_FSR) && ~isnan(fwhm_nm) && fwhm_nm > 0
        finesse = avg_FSR / fwhm_nm;
    else
        finesse = NaN;
    end

    all_results(i).Width_nm   = actual_widths_nm(i);
    all_results(i).Q_factor   = Q_factor;
    all_results(i).ER_Drop_dB = ER_drop_db;
    all_results(i).Finesse    = finesse;
end

W  = [all_results.Width_nm];
Q  = [all_results.Q_factor];
ER = [all_results.ER_Drop_dB];
F  = [all_results.Finesse];

%% ——— Color ———
C = lines(7);
skyblue = C(1,:);   % [0 0.4470 0.7410]

%% ——— Figure ———
fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 24 8]);

% ========================= (a) Q-factor =========================
ax1 = subplot(1,3,1);
plot(W, Q, 'o-', 'Color', skyblue, 'MarkerFaceColor', skyblue, ...
    'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Ring width (nm)');
ylabel('Q-factor');
ylim([0.85*min(Q) 1.12*max(Q)]);
xlim([min(W)-30 max(W)+30]);

text(0.04, 0.96, '(a)', 'Units','normalized', ...
    'FontSize', 21, 'FontWeight','bold', ...
    'FontName','Times New Roman', 'VerticalAlignment','top');

set(ax1, 'FontName','Times New Roman', 'FontSize', 21, ...
    'TickDir','in', 'LineWidth', 0.9, ...
    'XMinorTick','on', 'YMinorTick','on', 'Box','on');

% ========================= (b) Drop ER =========================
ax2 = subplot(1,3,2);
plot(W, ER, 'o-', 'Color', skyblue, 'MarkerFaceColor', skyblue, ...
    'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Ring width (nm)');
ylabel('Drop extinction ratio (dB)');
ylim([0.85*min(ER) 1.12*max(ER)]);
xlim([min(W)-30 max(W)+30]);

text(0.04, 0.96, '(b)', 'Units','normalized', ...
    'FontSize', 21, 'FontWeight','bold', ...
    'FontName','Times New Roman', 'VerticalAlignment','top');

set(ax2, 'FontName','Times New Roman', 'FontSize', 21, ...
    'TickDir','in', 'LineWidth', 0.9, ...
    'XMinorTick','on', 'YMinorTick','on', 'Box','on');

% ========================= (c) Finesse =========================
ax3 = subplot(1,3,3);
plot(W, F, 'o-', 'Color', skyblue, 'MarkerFaceColor', skyblue, ...
    'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Ring width (nm)');
ylabel('Finesse');
ylim([0.85*min(F) 1.12*max(F)]);
xlim([min(W)-30 max(W)+30]);

text(0.04, 0.96, '(c)', 'Units','normalized', ...
    'FontSize', 21, 'FontWeight','bold', ...
    'FontName','Times New Roman', 'VerticalAlignment','top');

set(ax3, 'FontName','Times New Roman', 'FontSize', 21, ...
    'TickDir','in', 'LineWidth', 0.9, ...
    'XMinorTick','on', 'YMinorTick','on', 'Box','on');

% ——— Export ———
exportgraphics(fig, 'Fig3_width_parameters.png', 'Resolution', 600);
exportgraphics(fig, 'Fig3_width_parameters.pdf', 'ContentType', 'vector');
fprintf('Done — Fig3_width_parameters.png / .pdf\n');

%% ——— Helper ———
function fwhm = calculate_fwhm(x, y, half_max)
    cross_idx = find(diff(sign(y - half_max)) ~= 0);
    if length(cross_idx) >= 2
        x1 = interp1(y(cross_idx(1):cross_idx(1)+1), ...
                     x(cross_idx(1):cross_idx(1)+1), half_max);
        x2 = interp1(y(cross_idx(2):cross_idx(2)+1), ...
                     x(cross_idx(2):cross_idx(2)+1), half_max);
        fwhm = abs(x2 - x1);
    else
        fwhm = NaN;
    end
    if fwhm < 0.001 || fwhm > 100, fwhm = NaN; end
end