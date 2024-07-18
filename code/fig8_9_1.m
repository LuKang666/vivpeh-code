%% fig.8_9 (位移，速度，电压和功率的期望值，填充瀑布图)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig89-1.mat');
swEPSfigure;
sigma = [0, 0.005, 0.01];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); % E[q1^2]
for i = 1 : 3
    xx = sigma(i) * ones(1, length(delta) + 2);
    z1 = zeros(1, length(delta) + 2);
    z1(2 : end - 1) = Eq1(:, i)';
    yy = zeros(1, length(delta) + 2);
    yy(1) = delta(1);
    yy(end) = delta(end);
    yy(2 : end - 1) = delta;
    fill3(xx, yy, z1, 'b', 'FaceAlpha', 0.5);
    hold on;
end
hold off;
grid on;
xlabel('$\sigma$'); ylabel('$\delta$'); zlabel('$E[q_1^2]$'); xticks([0, 0.005, 0.01]);
set(gca, 'LineWidth', 3.5);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); % E[p1^2]
for i = 1:3
    xx = sigma(i) * ones(1, length(delta) + 2);
    z1 = zeros(1, length(delta) + 2);
    z1(2:end - 1) = Ep1(:, i)';
    yy = zeros(1, length(delta) + 2);
    yy(1) = delta(1);
    yy(end) = delta(end);
    yy(2:end - 1) = delta;
    fill3(xx, yy, z1, 'b', 'FaceAlpha', 0.5);
    hold on;
end

hold off;
grid on;
xlabel('$\sigma$'); ylabel('$\delta$'); zlabel('$E[p_1^2]$'); xticks([0, 0.005, 0.01]);
set(gca, 'LineWidth', 3.5);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); % E[V^2]

for i = 1:3
    xx = sigma(i) * ones(1, length(delta) + 2);
    z1 = zeros(1, length(delta) + 2);
    z1(2:end - 1) = EV(:, i)';
    yy = zeros(1, length(delta) + 2);
    yy(1) = delta(1);
    yy(end) = delta(end);
    yy(2:end - 1) = delta;
    fill3(xx, yy, z1, 'b', 'FaceAlpha', 0.5);
    hold on;
end

hold off;
grid on;
xlabel('$\sigma$'); ylabel('$\delta$'); zlabel('$E[\overline{V}^2]$'); xticks([0, 0.005, 0.01]);
set(gca, 'LineWidth', 3.5);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); % E[P^2]

for i = 1:3
    xx = sigma(i) * ones(1, length(delta) + 2);
    z1 = zeros(1, length(delta) + 2);
    z1(2:end - 1) = EP(:, i)';
    yy = zeros(1, length(delta) + 2);
    yy(1) = delta(1);
    yy(end) = delta(end);
    yy(2:end - 1) = delta;
    fill3(xx, yy, z1, 'b', 'FaceAlpha', 0.5);
    hold on;
end

hold off;
grid on;
xlabel('$\sigma$'); ylabel('$\delta$'); zlabel('$E[\overline{P}]$'); xticks([0, 0.005, 0.01]);
set(gca, 'LineWidth', 3.5);
set(gcf, 'position', [0, 0, 500, 510]);
