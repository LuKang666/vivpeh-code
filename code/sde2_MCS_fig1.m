clc; clear;
%% 基本参数
alpha = 0.045; gamma = 0.6667; xi = 0.0043; St = 0.21; C_L0 = 0.3;
% Dime = 0.06; hp = 0.0005; hs = 0.00025; b = 0.04; len = 0.20; L = 0.15; % 几何参数
Dime = 0.06; hp = 0.267e-3; hs = 0.635e-3; b = 32.5e-3; len = 267e-3; L = 203e-3;

% Ep = 62e9; Es = 45e9; pa = 1.205; pp = 7600; ps = 1780; pc = 1780; k3 = 1.0e5; % 材料参数(PZT)
Ep = 66e9; Es = 70e9; pa = 1.205; pp = 7800; ps = 2730;
d31 = -190e-12; e33 = 13.28e-9;
% d31 = 340e-12; e33 = 1300 * 8.854e-12; % 压电参数(PZT)
% Ep = 2.5e9; Es = 4.5e9; pa = 1.205; pp = 1780; ps = 1820; pc = 52; k3 = 1.0e5; % 材料参数(PVDF)
% d31 = 33e-12; e33 = 18 * 8.854e-12; % 压电参数(PVDF)

A = Es / Ep; B = hs / hp;
EI = Ep * b * hp ^ 3 * (A * B ^ 3 + 6 * B ^ 2 + 12 * B + 8) / 12;
k_eq = 3 * EI / len ^ 3;
pb = (2 * hp * pp + hs * ps) / (2 * hp + hs);
m = b * pb * (2 * hp + hs);
m_eq = 33 * m * len / 140;
omn = sqrt(k_eq / m_eq);
Cp = 2 * e33 * b * len / hp;
% Cp = 2 * 20.3e-9;
R_eq = 2.46 * 1.0e6; % 负载阻抗
% lambda = 0.2631 * len / Dime;
% beta_nondim = 1 / R_eq / Cp / omn;
% beta = beta_nondim * omn;


%% 无量纲参数
lambda = 0.1;
beta = 0.4709;
mu = 0.0005;
% eta = 0;
% par = 9.0584;

% mu = pa * b * Dime ^ 2/8 / pi ^ 2 / m_eq / St ^ 2;
eta = Cp / m_eq / d31^2 / omn ^ 2;
delta = 0.0;
zeta = 2 * xi * omn + omn ^ 2 * beta / eta / (beta ^ 2 + omn ^ 2);
om1 = omn * sqrt(1 + omn ^ 2 / eta / (beta ^ 2 + omn ^ 2));
oms = (1 + delta) * om1; % 脱落频率
U0 = oms * Dime / 2 / pi / St; % 平均风速
par = 140 * U0 * pa / 33 / len / pb / (2*hp+hs);
% par = 140 * U0 / 33 / len / (2 * hp + hs);

sigma = 0.001; % 噪声强度
% sigma = 0.05;
% sigma = 5.7027e-11;
% sigma = 0.0;


%% 设置初值和时间序列
% x0 = 0.001 * ones(1, 4);
% x0 = 0.01 * ones(1, 4);
x0 = 0.1 * ones(1, 4);
% x0 = 0.5 * ones(1, 4);
dt = 0.0002;
% dt = 0.01;

sysn = @(t, x, w)[x(3); x(4);
            -om1 * om1 * x(1) - zeta * x(3) + mu * oms * oms * x(2) + par * x(2) * w(1);
            -oms * oms * x(2) + alpha * oms * x(4) - gamma / oms * x(4) * x(4) * x(4) + lambda * omn * x(3)]; %原受噪声的运动系统 Eq.(12)

timb = 200000; timt = 600000; %timb前的时间序列作为瞬态，不收集统计
tge = 100; hz = timb:tge:timt;
lenhz = length(hz); %对时间序列，每隔tge个间隔收集统计
Nsamp = 10; %样本数

Tim = 0;
%% 蒙特卡罗模拟
for k = 1:Nsamp
    disp([num2str(k), ' / ', num2str(Nsamp)]);
    noise = whitenoise(dt, timt + 1, sigma, 'WhiteNoise');
    % noise_Const_White = noise;
    % noise_Const_Daven = noise;
    noise_Fluct_White = noise;
    % noise_Fluct_Daven = noise;
    [t, x] = systemresponse(sysn, noise, dt, x0);
    Tim = t;
    % x_Const_White = x;
    % x_Const_Daven = x;
    x_Fluct_White = x;
    % x_Fluct_White_5 = x;
    % x_Fluct_Daven = x;
    % figure
    % subplot(1, 2, 1); plot(t, x(:, 1));
    % subplot(1, 2, 2); plot(t, x(:, 2));
    Y(:, 1) = 0.5 * (x(:, 3) .^ 2 + om1 * om1 * x(:, 1) .^ 2); %H1
    Y(:, 2) = 0.5 * (x(:, 4) .^ 2 + oms * oms * x(:, 2) .^ 2); %H2
    Y(:, 3) = mod(anglexy(om1 * x(:, 1), -x(:, 3)) - anglexy(oms * x(:, 2), -x(:, 4)) + pi, 2 * pi) - pi; % 相角差
    be = (k - 1) * lenhz + 1; ed = k * lenhz;
    YY(be:ed, 1:4) = x(hz, :); %稳态运动状态
    YY(be:ed, 5:7) = Y(hz, :); %稳态H1,H2,相角差
end

H1H2Fi = YY(:, 5:7); %提取出H1,H2,Psi数据
XV = YY(:, 1:4); %提取出稳态响应

%% 绘制结果图
npoints = 30;
figure;
[H1, H2, pH12] = getpdf2_1(YY(:, [5, 6]), [0, 1], [170, 190], [npoints, npoints]);
subplot(1, 3, 1); surf(H1, H2, pH12');
[H1, Fi, pH1F] = getpdf2_1(YY(:, [5, 7]), [0, 1], [-1.4, -0.8], [npoints, npoints]);
subplot(1, 3, 2); surf(H1, Fi, pH1F');
[H2, Fi, pH2F] = getpdf2_1(YY(:, [6, 7]), [170, 190], [-1.4, -0.8], [npoints, npoints]);
subplot(1, 3, 3); surf(H2, Fi, pH2F');

figure;
[H_1, pH1] = getpdf_1(YY(:, 5), [0, 1], npoints);
subplot(1, 3, 1); plot(H_1, pH1, 'ok-');
[H_2, pH2] = getpdf_1(YY(:, 6), [170, 190], npoints);
subplot(1, 3, 2); plot(H_2, pH2, 'ok-');
[F_i, pFi] = getpdf_1(YY(:, 7), [-1.4, -0.8], npoints);
subplot(1, 3, 3); plot(F_i, pFi, 'ok-');


fileDir = 'D:/Program Files/MATLAB/work/Script/VIVPEHs/data/';
savePath_1 = strcat(fileDir, 'Fig', num2str(67), '.mat');
savePath_2 = strcat(fileDir, 'Fig', num2str(67), '_nonres.mat'); % 对应非共振情况
savePath_3 = strcat(fileDir, 'Fig', num2str(3), '.mat');
savePath_4 = strcat(fileDir, 'Fig', num2str(67), '_nonres_whitenoise.mat');
save(savePath_1, '-append')
% save(savePath_4)
% save(savePath_2, 'H1', 'H2', 'Fi', 'H_1', 'H_2', 'F_i', 'pH12', 'pH1F', 'pH2F', 'pH1', 'pH2', 'pFi', 'H1H2Fi', 'XV', 'omn', 'om1', 'oms', '-append');
% save(savePath_3, 'Tim', 'timb', 'x_Fluct_White');

% [q1s, pq1s] = getpdf_1(XV(:, 1));
% [p1s, pp1s] = getpdf_1(XV(:, 3));
% [~, ~, pq1p1s] = getpdf2_1(XV(:, [1, 3]));
% intFun_q1 = q1s .^ 2 .* pq1s;
% intFun_p1 = p1s .^ 2 .* pp1s;
% intFun_q1p1 = q1s .* p1s .* pq1p1s;
% Eq1 = trapz(q1s, intFun_q1);
% Ep1 = trapz(p1s, intFun_p1);
% Eq1p1 = trapz(p1s, trapz(q1s, intFun_q1p1), 2);
% EV = (omn^2 / eta / (beta^2 + omn^2))^2 * Eq1 + (beta / eta / (beta^2 + omn^2))^2 * Ep1 + (2 * beta * omn^2 / eta^2 / (beta^2 + omn^2)^2) * Eq1p1;
% EP = eta * beta_nondim * EV;

% fileDir = 'D:/Program Files/MATLAB/work/Script/VIVPEHs/data/';
% savePath = strcat(fileDir, 'Fig', num2str(89), '.mat');
% save(savePath, 'Eq1', 'Ep1', 'EV', 'EP');
