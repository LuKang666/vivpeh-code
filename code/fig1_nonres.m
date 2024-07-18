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
% beta = 1 / R_eq / Cp;

%% 无量纲参数
lambda = 0.1;
beta = 0.4709;
mu = 0.0005;
% eta = 0;

% mu = pa * b * Dime ^ 2 / 8 / pi ^ 2 / m_eq / St ^ 2;
eta = m_eq * d31 ^ 2 * omn ^ 2 / Cp;
delta = 0.5;
zeta = 2 * xi * omn + omn ^ 2 * eta * beta / (beta ^ 2 + omn ^ 2);
om1 = omn * sqrt(1 + eta * omn ^ 2 / (beta ^ 2 + omn ^ 2));
oms = (1 + delta) * om1; % 脱落频率
U0 = oms * Dime / 2 / pi / St; % 平均风速

sigma = 0.001; % 噪声强度
% sigma = 0.005;
% sigma = 5.7027e-11;

r = 0.005; % 地面粗糙度
v10 = 20; % 参考风速
bar_w = @(w) 600*w/pi/v10;
S_w = @(w) 8*pi*r*v10^2*bar_w(w)^2 / (w*sign(1+bar_w(w)^2)*abs(1+bar_w(w)^2)^(4/3)+eps);
% S_1 = sqrt(2 * sigma) * S_w(om1 - oms); S_2 = sqrt(2 * sigma) * S_w(om1 + oms);
S_1 = sigma / 2 / pi; S_2 = sigma / 2 / pi;

c = 140 * pa * U0 / 33 / pb / len / (2 * hp + hs);
a1 = @(at) -zeta * at(1) + pi * c^2 * at(2) * (S_1 + S_2) / 2 / oms / oms;
a2 = @(at) alpha * oms * at(2) - 3 * gamma * at(2)^2 / 2 / oms;
b11 = @(at) pi * c^2 * at(1) * at(2) * (S_1 + S_2) / oms / oms;

drift = {a1; a2}; %漂移系数矩阵
diffu = {b11, "0"; "0", "0"}; %扩散系数矩阵
% domain = [0.001, 3; 170, 220; -1.4, 0]; %离散域
domain_H = [0.001, 1.0; 170, 190];
% domain_H = [0.001, 2.0; 1, 380; -pi, pi];
% domain_H = [0.0038, 1.5016; 1.0154, 367.1; -pi, pi];
class = {'zero', 'zero'; 'zero', 'zero'}; %free是自由边界,period是周期边界,zero是零值边界
NH = [100, 100]; %输入求解域各维划分份数，两边界上的点也编号，比如分成10份，那么10-1共9个区间
% NH = [60, 23, 54];
% NH = [64, 36, 54]; % 对应Fig67.mat
% NH = [14, 14, 196]; % 对应Fig67_1.mat
% p = getPH(drift, diffu, domain, N, class);
% p = getPH2(drift, diffu, domain, N, class);
[A, C, v] = getA(drift, diffu, domain_H, NH, class);
NH = NH(:)'; Dim = length(NH); dx = diff(domain_H') ./ (NH - 1); dV = prod(dx); Points = prod(NH);

% 以下求瞬态
AA = A(v, v);
sys = @(t, p)AA * p(:);
p0 = zeros(Points, 1); p0(v) = 1; p0 = p0 / (C * p0 * dV); p0 = p0(v);
[t, pv] = ode45(sys, [0, 5], p0); % Transient state at enough time
figure;
plot(t, pv(:, end));
p = zeros(Points, 1); p(v) = pv(end, :);
p = p / (C * p * dV); %probability normalization
p = abs(p);
% showpxpxy(p, N, domain, C);

[H, pH, pHH] = showpxpxy(p, NH, domain_H, C, [1, 2], [1, 2]);
figure;
surf(H{1}, H{2}, pHH{1}');

figure;
subplot(1, 2, 1); plot(H{1}, pH{1}, 'k-');
subplot(1, 2, 2); plot(H{2}, pH{2}, 'k-');

fileDir = 'D:\Program Files\MATLAB\work\Script\VIVPEHs\data\';
% savePath = strcat(fileDir, 'Fig', num2str(67), '.mat'); % 对应共振情况
savePath_2 = strcat(fileDir, 'Fig', num2str(67), '_nonres.mat'); % 对应非共振情况
save(savePath_2, 'H', 'pH', 'pHH', 'drift', 'diffu', 'domain_H', 'class', 'NH', 'p', 'C');
