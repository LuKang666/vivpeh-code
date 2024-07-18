%% fig.8_9 (位移，速度，电压和功率的期望值)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig89-1.mat');
swEPSfigure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % E[q1^2] & E[p1^2]
S(1) = subplot(1, 2, 1);
plot(delta', Eq1(:, 1), '-or'); hold on;
plot(delta', Eq1(:, 2), '-^b'); hold on;
plot(delta', Eq1(:, 3), '-sg');
ylim([0, 2.0e-5]);
xlabel('$\delta$'); ylabel('$E[q_1^2]$');
set(legend('$\sigma=0$', ...
    '$\sigma=0.005$', ...
    '$\sigma=0.01$', ...
    'Location', 'northwest', 'FontSize', 14));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');
S(2) = subplot(1, 2, 2);
plot(delta', Ep1(:, 1), '-or'); hold on;
plot(delta', Ep1(:, 2), '-^b'); hold on;
plot(delta', Ep1(:, 3), '-sg');
ylim([0, 0.21])
xlabel('$\delta$'); ylabel('$E[p_1^2]$');
set(legend('$\sigma=0$', ...
    '$\sigma=0.005$', ...
    '$\sigma=0.01$', ...
    'Location', 'northwest', 'FontSize', 14));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % E[V^2] & E[P]
S(1) = subplot(1, 2, 1);
plot(delta', EV(:, 1), '-or'); hold on;
plot(delta', EV(:, 2), '-^b'); hold on;
plot(delta', EV(:, 3), '-sg');
ylim([0, 6.5e-5]);
xlabel('$\delta$'); ylabel('$E[\overline{V}^2]$');
set(legend('$\sigma=0$', ...
    '$\sigma=0.005$', ...
    '$\sigma=0.01$', ...
    'Location', 'northwest', 'FontSize', 14));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');
S(2) = subplot(1, 2, 2);
plot(delta', EP(:, 1), '-or'); hold on;
plot(delta', EP(:, 2), '-^b'); hold on;
plot(delta', EP(:, 3), '-sg');
ylim([0, 3.0e-7]);
xlabel('$\delta$'); ylabel('$E[\overline{P}]$');
set(legend('$\sigma=0$', ...
    '$\sigma=0.005$', ...
    '$\sigma=0.01$', ...
    'Location', 'northwest', 'FontSize', 14));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');
