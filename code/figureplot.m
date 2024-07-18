%% fig.1 (能量边缘概率密度解析解与数值解) (数据保存在data文件夹中，运行程序前请先定位到该文件夹)
clc; clear;
load analysis-2.mat; load simulation-3.mat;
swEPSfigure;
figure;
subplot(1, 3, 1);
plot(H{1}, pH{1}, 'k-');
hold on;
plot(H_1, pH1, 'ok');
xlabel('$H_1$'); ylabel('$p(H_1)$');
set(legend('Analytical solution', ...
    'Simulation', ...
    'Location', 'best', 'FontSize', 12, 'LineWidth', 1));
set(gca, 'LineWidth', 3.5);

subplot(1, 3, 2);
plot(H{2}, pH{2}, 'k-', 'LineWidth', 2);
hold on;
plot(H_2, pH2, 'ok', 'LineWidth', 2, 'MarkerSize', 12);
ylim([0, 0.27]);
xlabel('$H_2$'); ylabel('$p(H_2)$');
set(legend('Analytical solution', ...
    'Simulation', ...
    'Location', 'best', 'FontSize', 12, 'LineWidth', 1));
set(gca, 'LineWidth', 3.5);

subplot(1, 3, 3);
plot(H{3}, pH{3}, 'k-', 'LineWidth', 2);
hold on;
plot(F_i, pFi, 'ok', 'LineWidth', 2, 'MarkerSize', 12);
ylim([0, 6.8]);
xlabel('$\psi$'); ylabel('$p(\psi)$');
set(legend('Analytical solution', ...
    'Simulation', ...
    'Location', 'best', 'FontSize', 12, 'LineWidth', 1));
set(gca, 'LineWidth', 3.5);



%% fig.2 (能量联合概率密度解析解与数值解)
clc; clear;
load analysis-2.mat; load simulation-3.mat;
swEPSfigure;
factor2 = max(max(pH1F)) / max(max(pHH{2}));
factor3 = max(max(pHH{3})) / max(max(pH2F));
pHH{2} = pHH{2} * factor2;
pH2F = pH2F * factor3;
figure;
subplot(2, 3, 1); surf(H{1}, H{2}, pHH{1}');
colormap('jet');
shading interp;
zlim([0, 1]);
set(gca, 'LineWidth', 3.5);
xlabel('$H_1$'); ylabel('$H_2$'); zlabel('$p(H_1,H_2)$');
view(45, 30);
subplot(2, 3, 2); surf(H{1}, H{3}, pHH{2}');
colormap('jet');
shading interp;
zlim([0, 30]);
ylim([-1.4, -0.8]);
yticks([-1.4, -1.1, -0.8]);
yticklabels({'-1.4', '-1.1', '-0.8'});
set(gca, 'LineWidth', 3.5);
xlabel('$H_1$'); ylabel('$\psi$'); zlabel('$p(H_1,\psi)$');
view(45, 30);
subplot(2, 3, 3); surf(H{2}, H{3}, pHH{3}');
colormap('jet');
shading interp;
zlim([0, 2]);
ylim([-1.4, -0.8]);
yticks([-1.4, -1.1, -0.8]);
yticklabels({'-1.4', '-1.1', '-0.8'});
set(gca, 'LineWidth', 3.5);
xlabel('$H_2$'); ylabel('$\psi$'); zlabel('$p(H_2,\psi)$');
view(45, 30);

subplot(2, 3, 4); surf(H1, H2, pH12');
colormap('jet');
zlim([0, 1]);
set(gca, 'LineWidth', 3.5);
xlabel('$H_1$'); ylabel('$H_2$'); zlabel('$p(H_1,H_2)$');
view(45, 30);
subplot(2, 3, 5); surf(H1, Fi, pH1F');
colormap('jet');
zlim([0, 30]);
ylim([-1.4, -0.8]);
yticks([-1.4, -1.1, -0.8]);
yticklabels({'-1.4', '-1.1', '-0.8'});
set(gca, 'LineWidth', 3.5);
xlabel('$H_1$'); ylabel('$\psi$'); zlabel('$p(H_1,\psi)$');
view(45, 30);
subplot(2, 3, 6); surf(H2, Fi, pH2F');
colormap('jet');
zlim([0, 2]);
ylim([-1.4, -0.8]);
yticks([-1.4, -1.1, -0.8]);
yticklabels({'-1.4', '-1.1', '-0.8'});
set(gca, 'LineWidth', 3.5);
xlabel('$H_2$'); ylabel('$\psi$'); zlabel('$p(H_2,\psi)$');
view(45, 30);



%% fig.3 (q1和p1的稳态时程图)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig3.mat');
swEPSfigure;
alpha1 = 0.1; alpha2 = 0.5; len = 100000; aa = 500; linewid_2 = 2;
figure(1);
h1 = plot3(Tim(timb:end), x_Fluct_White(timb:end, 1), x_Fluct_White(timb:end, 3));
% h1 = plot3(Tim(end-len:aa:end), x_Fluct_White(end-len:aa:end, 1), x_Fluct_White(end-len:aa:end, 3));
h1.Color(4) = alpha1;
set(gca, 'LineWidth', 3.5);
ylim([-0.04, 0.04]); zlim([-2, 2]);
xlabel('$t$'); ylabel('$q_1$'); zlabel('$p_1$');
hold on;
h3 = plot3(Tim(timb:end), x_Fluct_White(timb:end, 1), -2 * ones(size(x_Fluct_White(timb:end, 3))), 'r-'); % 绘制投影
% h3 = plot3(Tim(end - len:aa:end), x_Fluct_White(end - len:aa:end, 1), -2 * ones(size(x_Fluct_White(end - len:aa:end, 3))), 'r-', 'LineWidth', linewid_2); % 绘制投影
hold on;
h4 = plot3(Tim(timb:end), 0.04 * ones(size(x_Fluct_White(timb:end, 1))), x_Fluct_White(timb:end, 3), 'g-');
% h4 = plot3(Tim(end - len:aa:end), 0.04 * ones(size(x_Fluct_White(end - len:aa:end, 1))), x_Fluct_White(end - len:aa:end, 3), 'g-', 'LineWidth', linewid_2);
h3.Color(4) = alpha2;
h4.Color(4) = alpha2;
set(gca, 'Clipping', 'on', 'LooseInset', [0 0 0 0]);
set(gcf, 'Position', [0, 0, 500, 450]);

figure(2);
h2 = plot3(Tim(timb:end), x_Fluct_White(timb:end, 2), x_Fluct_White(timb:end, 4));
% h2 = plot3(Tim(end-len:aa:end), x_Fluct_White(end-len:aa:end, 2), x_Fluct_White(end-len:aa:end, 4));
h2.Color(4) = alpha1;
set(gca, 'LineWidth', 3.5);
ylim([-0.8, 0.8]); zlim([-40, 40]);
xlabel('$t$'); ylabel('$q_2$'); zlabel('$p_2$');
hold on;
h5 = plot3(Tim(timb:end), x_Fluct_White(timb:end, 2), -40 * ones(size(x_Fluct_White(timb:end, 4))), 'r-'); % 绘制投影
% h5 = plot3(Tim(end - len:aa:end), x_Fluct_White(end - len:aa:end, 2), -40 * ones(size(x_Fluct_White(end - len:aa:end, 4))), 'r-', 'LineWidth', linewid_2); % 绘制投影
hold on;
h6 = plot3(Tim(timb:end), 0.8 * ones(size(x_Fluct_White(timb:end, 2))), x_Fluct_White(timb:end, 4), 'g-');
% h6 = plot3(Tim(end - len:aa:end), 0.8 * ones(size(x_Fluct_White(end - len:aa:end, 2))), x_Fluct_White(end - len:aa:end, 4), 'g-', 'LineWidth', linewid_2);
h5.Color(4) = alpha2;
h6.Color(4) = alpha2;
set(gca, 'Clipping', 'on', 'LooseInset', [0 0 0 0]);
set(gcf, 'Position', [0, 0, 500, 450]);

%% fig.4 (共振时用高斯白噪声替代Davenport PSD)
clc; clear;
r = 0.005; % 地面粗糙度
v10 = 20; % 参考风速
freqs = logspace(-3, 0, 1000); % 频率向量
omega = 2 * pi * freqs;
om_bar = 600 * omega / v10 / pi;
S_omega = 8 * pi * r * v10 ^ 2 * om_bar .^ 2 ./ (omega .* sign(1 + om_bar .^ 2) .* abs(1 + om_bar .^ 2) .^ (4/3) + eps); % 目标功率谱密度（加了eps以避免分母为零）
f_bar = 1200 * freqs / v10;
S_freq = 4 * r * v10 ^ 2 * f_bar .^ 2 ./ (freqs .* sign(1 + f_bar .^ 2) .* abs(1 + f_bar .^ 2) .^ (4/3) + eps); % 目标功率谱密度（加了eps以避免分母为零）
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
figure;
loglog(omega, S_omega, 'k-', 'LineWidth', 2);
[sigma, index] = max(S_omega);
hold on;
loglog(omega, sigma * ones(1, length(omega)), 'k--', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20, 'LineWidth', 3.5, 'FontWeight', 'bold');
xlabel('$\omega$'); ylabel('$S(\omega)$');
ylim([0,500]);
text(omega(1), sigma, '$$\frac{\sigma}{2\pi}$$', 'interpreter', 'latex', 'fontsize', 20)
legend('Davenport PSD', 'WhiteNoise', 'Location', 'northeast');



%% fig.5 (X1的功率谱密度)
clc; clear;
load simulation-4.mat;
swEPSfigure;
figure;
% 使用MATLAB的fft函数进行傅立叶变换
N = length(x_Const_White(:, 1)); % 数据点数
X1 = fft(x_Const_White(:, 1));
% 计算单边频谱
X1_single = X1(1:N / 2 + 1);
% 计算功率谱密度
psd1 = (1 / (N ^ 2)) * abs(X1_single) .^ 2;
% 构建频率向量
fs = 1000; % 假设采样频率为1000 Hz
f = (0:length(X1_single) - 1) * (fs / N);
omega = 2 * pi * f;
% 绘制功率谱密度随角频率的变化图像
plot(omega, psd1, 'k-')
set(gca, 'LineWidth', 3.5);
xlim([12, 12.5]);
xlabel('Frequency $$\omega$$ [rad/s]')
ylabel('PSD of displacement $$X_1$$')
hold on
X2 = fft(x_Fluct_White_5(:, 1));
X2_single = X2(1:N / 2 + 1);
psd2 = (1 / (N ^ 2)) * abs(X2_single) .^ 2;
plot(omega, psd2, 'k--')
set(legend('$\sigma=0$', '$\sigma=0.005$', 'Location', 'best', 'FontSize', 12, 'LineWidth', 1));
set(gca, 'LooseInset', [0 0 0 0]);
set(gcf, 'Position', [50, 100, 450, 400]);




