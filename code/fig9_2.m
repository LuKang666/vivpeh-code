%% fig.9_2 (共振时电压和功率的矩随eta的变化,sigma=0.01, beta_nondim=0.008)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig89-3.mat');
% load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig89-4.mat');
swEPSfigure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
span = 5; poly = 1;
% EV_smooth = smooth(eta, EV, span, 'sgolay', poly);
EP_smooth = smooth(eta, EP, span, 'sgolay', poly);
% EV_smooth = smooth(eta, EV, span, 'rloess');
% EP_smooth = smooth(eta, EP, span, 'rloess');

yyaxis left;
plot(eta, EV, '-ob');
xlabel('$\eta$'); ylabel('$E[\overline{V}^2]$');
yyaxis right;
plot(eta, EP_smooth, '-sr');
xlim([0, 3.1]);
xlabel('$\eta$'); ylabel('$E[\overline{P}]$');
set(legend('$E[\overline{V}^2]$', ...
    '$E[\overline{P}]$', ...
    'Location', 'northeast', 'FontSize', 14));
% set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.4, 0.2]);
set(gcf, 'position', [50, 100, 650, 500]);
set(gca, 'LineWidth', 3.5);
