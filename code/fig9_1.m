%% fig.9_1 (共振时电压和功率的矩随beta_nondim的变化, sigma=0.01, eta=0.526)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig89-2.mat');
swEPSfigure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
span = 5;
EP_smooth = smooth(beta_nondim, EP, span, 'sgolay', 1);
yyaxis left;
plot(beta_nondim, EV, '-ob');
xlabel('$\overline{\beta}$'); ylabel('$E[\overline{V}^2]$');
yyaxis right;
plot(beta_nondim, EP_smooth, '-sr');
xlim([0, 3.1]);
xlabel('$\overline{\beta}$'); ylabel('$E[\overline{P}]$');
set(legend('$E[\overline{V}^2]$', ...
    '$E[\overline{P}]$', ...
    'Location', 'northeast', 'FontSize', 14));
% set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.4, 0.2]);
set(gcf, 'position', [50, 100, 650, 500]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');
