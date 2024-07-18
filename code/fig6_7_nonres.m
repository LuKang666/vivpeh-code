%% fig.6_7 (非共振时响应边缘概率密度和联合概率密度的解析解与数值解)
clc; clear;
data = load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig67_nonres.mat'); % 对应非共振情况 delta = 0.5
data2 = load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig67_nonres_whitenoise.mat'); % 非共振高斯白噪声
swEPSfigure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.001;
zeta = data.zeta; oms = data.oms; om1 = data.om1; alpha = data.alpha; gamma = data.gamma;
c = data.par * sqrt(2 * sigma);
r = 0.005; % 地面粗糙度
v10 = 20; % 参考风速
bar_w = @(w) 600 * w / pi / v10;
S_w = @(w) 8 * pi * r * v10 ^ 2 * bar_w(w) ^ 2 / (w * sign(1 + bar_w(w) ^ 2) * abs(1 + bar_w(w) ^ 2) ^ (4/3) + eps);
S_1 = sqrt(2 * sigma) * S_w(om1 - oms); S_2 = sqrt(2 * sigma) * S_w(om1 + oms);
% S_1 = sigma / 2 / pi; S_2 = sigma / 2 / pi;

Num = 60;
q1 = linspace(min(data.XV(:, 1)), max(data.XV(:, 1)), Num); p1 = linspace(min(data.XV(:, 3)), max(data.XV(:, 3)), Num);
q2 = linspace(min(data.XV(:, 2)), max(data.XV(:, 2)), Num); p2 = linspace(min(data.XV(:, 4)), max(data.XV(:, 4)), Num);
pq1p1 = zeros(Num); kappa = zeros(Num);
for i = 1:Num
    for j = 1:Num
        h1 = 0.5 * p1(j) * p1(j) + 0.5 * om1 * om1 * q1(i) * q1(i);
        h2 = 0.5 * p2(j) * p2(j) + 0.5 * oms * oms * q2(i) * q2(i);
        kappa(i, j) = log(2 * zeta * oms^2 * h1 - pi * c^2 * h2 * (S_1 + S_2)) + log(3 * gamma * h2^2 - 2 * alpha * oms^2 * h2);
        pq1p1(i, j) = (om1 / 2 / pi) * exp(-real(kappa(i, j)));
    end
end
dq1 = (max(q1) - min(q1)) / (Num - 1);
dp1 = (max(p1) - min(p1)) / (Num - 1);
pq1p1 = pq1p1 / sum(pq1p1 * dq1 * dp1, 'all'); % 归一化

[q1s, p1s, pq1p1s] = getpdf2(data.XV(:, [1, 3]), [], [], [Num, Num]);
[q1s_w, p1s_w, pq1p1s_w] = getpdf2(data2.XV(:, [1, 3]), [], [], [Num, Num]);
pq1p1s_w = pq1p1s_w * (max(max(pq1p1s)) / max(max(pq1p1s_w)));
constant1 = trapz(p1, trapz(q1, pq1p1), 2);
% constant2 = trapz(p1, trapz(q1, pq1p1s), 2);
pq1p1 = pq1p1 * (max(max(pq1p1s)) / max(max(pq1p1)));
figure(2); % 二维密度p(q1,p1)
% cv = gray; colormap(cv(floor(end * 0.4):end, :));
surf(q1s, p1s, pq1p1s');
colormap('cool');
hold on;
surf(q1, p1, -200 + 0 * pq1p1s', pq1p1s');
set(gca, 'zlim', [-200, 400], 'LineWidth', 3.5, ...
    'view', [-37.5, 15]);
xlabel('$q_1$'); ylabel('$p_1$'); zlabel('PDF $p(q_1,p_1)$'); zticks([0, 200, 400]);
set(gcf, 'position', [0, 0, 500, 510]);

% S(2) = subplot(1, 2, 1); surf(q1, p1, pq1p1');
figure(1);
surf(q1s_w, p1s_w, pq1p1s_w');
hold on;
surf(q1, p1, zeros(Num));
colormap('cool');
shading interp;
hold on;
surf(q1, p1, -200 + 0 * pq1p1s_w', pq1p1s_w');
shading interp;
set(gca, 'xlim', [min(data.XV(:, 1)), max(data.XV(:, 1))], ...
    'ylim', [min(data.XV(:, 3)), max(data.XV(:, 3))], ...
    'zlim', [-200, 400], 'LineWidth', 3.5, ...
    'view', [-37.5, 15]);
xlabel('$q_1$'); ylabel('$p_1$'); zlabel('PDF $p(q_1,p_1)$'); zticks([0, 200, 400]);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % 一维密度p(q1),p(p1)
[q1s, pq1s] = getpdf_1(data.XV(:, 1));
[q1s_w, pq1s_w] = getpdf_1(data2.XV(:, 1));
pq1s_w = pq1s_w * (max(pq1s) / max(pq1s_w));
pq1 = trapz(p1, pq1p1, 2); pq1 = reshape(pq1, [Num, 1]);
pq1 = pq1 * (max(pq1s) / max(pq1));
% S(1) = subplot(1, 2, 1); plot(q1s, pq1s, 'ok', q1, pq1, 'k');
figure(3);
plot(q1s_w, pq1s_w, 'k', q1s, pq1s, 'ok');
ylim([0, 135]);
xlabel('$q_1$'); ylabel('PDF $p(q_1)$');
legend('from analytical solution', 'from simulation', 'Location', 'northeast', 'FontSize', 14);
set(gcf, 'position', [0, 0, 500, 550]);
set(gca, 'LineWidth', 3.5);
[p1s, pp1s] = getpdf_1(data.XV(:, 3));
[p1s_w, pp1s_w] = getpdf_1(data2.XV(:, 3));
pp1s_w = pp1s_w * (max(pp1s) / max(pp1s_w));
pp1 = trapz(q1, pq1p1, 1); pp1 = reshape(pp1, [Num, 1]);
pp1 = pp1 * (max(pp1s) / max(pp1));
% S(2) = subplot(1, 2, 2); plot(p1s, pp1s, 'ok', p1, pp1, 'k');
figure(4); plot(p1s_w, pp1s_w, 'k', p1s, pp1s, 'ok');
ylim([0, 2.3]);
xlabel('$p_1$'); ylabel('PDF $p(p_1)$');
legend('from analytical solution', 'from simulation', 'Location', 'northeast', 'FontSize', 14);
set(gcf, 'position', [0, 0, 500, 550]);
set(gca, 'LineWidth', 3.5);
