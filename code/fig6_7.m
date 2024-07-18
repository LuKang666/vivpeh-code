%% fig.6_7 (共振时响应边缘概率密度和联合概率密度的解析解与数值解)
clc; clear;
data = load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig67.mat'); % 对应共振情况
% data = load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig67_nonres.mat'); % 对应非共振情况
swEPSfigure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_H1 = data.NH(1); Num_H2 = data.NH(2);
q1 = linspace(min(data.XV(:, 1)), max(data.XV(:, 1)), Num_H1); p1 = linspace(min(data.XV(:, 3)), max(data.XV(:, 3)), Num_H1);
q2 = linspace(min(data.XV(:, 2)), max(data.XV(:, 2)), Num_H2); p2 = linspace(min(data.XV(:, 4)), max(data.XV(:, 4)), Num_H2);
pq1p1 = zeros(Num_H1);

for i = 1:Num_H1

    for j = 1:Num_H1
        h1 = 0.5 * (p1(j) ^ 2 + data.om1 ^ 2 * q1(i) ^ 2);
        C = abs(data.H{1} - h1);
        C_min = min(C);
        index = find(C == C_min); % 找出H{1}中与h1最接近的元素的下标
        pq1p1(i, j) = data.om1 / 2 / pi * data.pH{1}(index);
    end

end

points = 60;
% points = 30;
% [q1s, p1s, pq1p1s] = getpdf2_1(data.XV(:, [1, 3]), [], [], [points, points]);
[~, ~, pq1p1s] = getpdf2(data.XV(:, [1, 3]), [], [], [points, points]);
pq1p1s = 2 * pq1p1s;
[m, n] = size(pq1p1s);
pq1p1s(1, :) = 0;
pq1p1s(m, :) = 0;
pq1p1s(:, 1) = 0;
pq1p1s(:, n) = 0;

figure(2); % 二维密度p(q1,p1)
% cv = gray; colormap(cv(floor(end * 0.4):end, :));
surf(q1, p1, pq1p1s');
colormap('cool');
hold on;
surf(q1, p1, -30 + 0 * pq1p1s', pq1p1s');
set(gca, 'xlim', [min(data.XV(:, 1)), max(data.XV(:, 1))], ...
         'ylim', [min(data.XV(:, 3)), max(data.XV(:, 3))], ...
         'zlim', [-30, 50], ...
         'LineWidth', 3.5, ...
         'view', [-37.5, 30], ...
         'LooseInset', [0 0 0 0]);
xlabel('$q_1$'); ylabel('$p_1$'); zlabel('PDF $p(q_1,p_1)$'); zticks([0, 25, 50]);
set(gcf, 'position', [0, 0, 500, 510]);

figure(1);
surf(q1, p1, pq1p1');
colormap('cool'); shading interp;
hold on;
surf(q1, p1, -30 + 0 * pq1p1', pq1p1');
shading interp;
set(gca, 'xlim', [min(data.XV(:, 1)), max(data.XV(:, 1))], ...
    'ylim', [min(data.XV(:, 3)), max(data.XV(:, 3))], ...
    'zlim', [-30, 50], ...
    'LineWidth', 3.5, ...
    'view', [-37.5, 30], ...
    'LooseInset', [0 0 0 0]);
xlabel('$q_1$'); ylabel('$p_1$'); zlabel('PDF $p(q_1,p_1)$'); zticks([0, 25, 50]);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); % 一维密度p(q1),p(p1)
smooth_factor = 100;
[q1s, pq1s] = getpdf_1(data.XV(:, 1));
pq1 = trapz(p1, pq1p1, 2); pq1 = reshape(pq1, [Num_H1, 1]);
% pq1_smooth = smooth(q1, pq1, smooth_factor, 'sgolay');
plot(q1, pq1, 'k', q1s, pq1s, 'ok');
ylim([0, 55]);
xlabel('$q_1$'); ylabel('PDF $p(q_1)$');
legend('from analytical solution', 'from simulation', 'Location', 'northeast', 'FontSize', 14);
set(gcf, 'position', [0, 0, 500, 550]);
set(gca, 'LineWidth', 3.5);
[p1s, pp1s] = getpdf_1(data.XV(:, 3));
pp1 = trapz(q1, pq1p1, 1); pp1 = reshape(pp1, [Num_H1, 1]);
% pp1_smooth = smooth(p1, pp1, smooth_factor, 'sgolay');
figure(4);
plot(p1, pp1, 'k', p1s, pp1s, 'ok');
ylim([0, 0.9]);
xlabel('$p_1$'); ylabel('PDF $p(p_1)$');
legend('from analytical solution', 'from simulation', 'Location', 'northeast', 'FontSize', 14);
set(gcf, 'position', [0, 0, 500, 550]);
set(gca, 'LineWidth', 3.5);
