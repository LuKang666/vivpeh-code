%% fig9_1 (电压和功率的矩随beta_nondim的变化)
clc; clear;
%% 基本参数
alpha = 0.045; gamma = 0.6667; xi = 0.0043; St = 0.21; C_L0 = 0.3;
Dime = 0.06; hp = 0.267e-3; hs = 0.635e-3; b = 32.5e-3; len = 267e-3; L = 203e-3;
Ep = 66e9; Es = 70e9; pa = 1.205; pp = 7800; ps = 2730;
d31 = -190e-12; e33 = 13.28e-9;

A = Es / Ep; B = hs / hp;
EI = Ep * b * hp ^ 3 * (A * B ^ 3 + 6 * B ^ 2 + 12 * B + 8) / 12;
k_eq = 3 * EI / len ^ 3;
pb = (2 * hp * pp + hs * ps) / (2 * hp + hs);
m = b * pb * (2 * hp + hs);
m_eq = 33 * m * len / 140;
omn = sqrt(k_eq / m_eq);
Cp = 2 * e33 * b * len / hp;
R_eq = 2.46 * 1.0e6; % 负载阻抗
% beta_nondim = 1 / R_eq / Cp / omn;
beta_nondim = linspace(0.1, 3, 20);
beta = beta_nondim * omn;

%% 无量纲参数
lambda = 0.1;
% beta = 0.4709;
mu = 0.0005;
eta = 0.526;
% par = 9.0584;

% eta = Cp / m_eq / d31 ^ 2 / omn ^ 2;
zeta = 2 * xi * omn + omn ^ 2 * beta ./ eta ./ (beta .^ 2 + omn ^ 2);
om1 = omn * sqrt(1 + omn ^ 2 / eta ./ (beta .^ 2 + omn ^ 2));

delta = 0.0;
oms = (1 + delta) .* om1; % 脱落频率
U0 = oms * Dime / 2 / pi / St; % 平均风速
par = 140 * U0 * pa / 33 / len / pb / (2 * hp + hs);

sigma = 0.01; % 噪声强度
% sigma = 0.0;

%% 设置初值和时间序列
x0 = 0.1 * ones(1, 4);
dt = 0.0002;

Eq1 = zeros(length(beta), 1);
Ep1 = zeros(length(beta), 1);
Eq1p1 = zeros(length(beta), 1);
EV = zeros(length(beta), 1);
EP = zeros(length(beta), 1);

for ii = 1:length(beta)
    disp(['ii: ', num2str(ii), ' / ', num2str(length(beta))]);
    sysn = @(t, x, w)[x(3); x(4);
                      -om1(ii) * om1(ii) * x(1) - zeta(ii) * x(3) + mu * oms(ii) * oms(ii) * x(2) + par(ii) * x(2) * w(1);
                      -oms(ii) * oms(ii) * x(2) + alpha * oms(ii) * x(4) - gamma / oms(ii) * x(4) * x(4) * x(4) + lambda * omn * x(3)]; %原受噪声的运动系统 Eq.(12)

    % timb = 200000; timt = 600000; %timb前的时间序列作为瞬态，不收集统计
    timb = 100000; timt = 200000;
    tge = 100; hz = timb:tge:timt;
    lenhz = length(hz); %对时间序列，每隔tge个间隔收集统计
    Nsamp = 10; %样本数

    %% 蒙特卡罗模拟
    for k = 1:Nsamp
        disp([num2str(k), ' / ', num2str(Nsamp)]);
        noise = whitenoise(dt, timt + 1, sigma, 'WhiteNoise');
        % noise_Const_White = noise;
        % noise_Const_Daven = noise;
        noise_Fluct_White = noise;
        % noise_Fluct_Daven = noise;
        [~, x] = systemresponse(sysn, noise, dt, x0);
        % x_Const_White = x;
        % x_Const_Daven = x;
        x_Fluct_White = x;
        % x_Fluct_White_5 = x;
        % x_Fluct_Daven = x;
        % figure
        % subplot(1, 2, 1); plot(t, x(:, 1));
        % subplot(1, 2, 2); plot(t, x(:, 2));
        Y(:, 1) = 0.5 * (x(:, 3) .^ 2 + om1(ii) * om1(ii) * x(:, 1) .^ 2); %H1
        Y(:, 2) = 0.5 * (x(:, 4) .^ 2 + oms(ii) * oms(ii) * x(:, 2) .^ 2); %H2
        Y(:, 3) = mod(anglexy(om1(ii) * x(:, 1), -x(:, 3)) - anglexy(oms(ii) * x(:, 2), -x(:, 4)) + pi, 2 * pi) - pi; % 相角差
        be = (k - 1) * lenhz + 1; ed = k * lenhz;
        YY(be:ed, 1:4) = x(hz, :); %稳态运动状态
        YY(be:ed, 5:7) = Y(hz, :); %稳态H1,H2,相角差
    end

    H1H2Fi = YY(:, 5:7); %提取出H1,H2,Psi数据
    XV = YY(:, 1:4); %提取出稳态响应

    [q1s, pq1s] = getpdf_1(XV(:, 1));
    [p1s, pp1s] = getpdf_1(XV(:, 3));
    [~, ~, pq1p1s] = getpdf2_1(XV(:, [1, 3]));
    intFun_q1 = q1s .^ 2 .* pq1s;
    intFun_p1 = p1s .^ 2 .* pp1s;
    intFun_q1p1 = q1s .* p1s .* pq1p1s;
    Eq1(ii) = trapz(q1s, intFun_q1);
    Ep1(ii) = trapz(p1s, intFun_p1);
    Eq1p1(ii) = trapz(p1s, trapz(q1s, intFun_q1p1), 2);
    EV(ii) = (omn ^ 2 / eta ./ (beta(ii) .^ 2 + omn ^ 2)) ^ 2 * Eq1(ii) + (beta(ii) / eta / (beta(ii) ^ 2 + omn ^ 2)) ^ 2 * Ep1(ii) + (2 * beta(ii) * omn ^ 2 / eta ^ 2 / (beta(ii) ^ 2 + omn ^ 2) ^ 2) * Eq1p1(ii);
    EP(ii) = eta * beta_nondim(ii) * EV(ii);

end

fileDir = 'D:/Program Files/MATLAB/work/Script/VIVPEHs/data/';
savePath = strcat(fileDir, 'Fig', num2str(89), '-', num2str(2), '.mat'); % 共振情形
savePath_2 = strcat(fileDir, 'Fig', num2str(89), '-', num2str(2), '_nonres.mat'); % 非共振情形 delta=0.5
% save(savePath_2, 'beta_nondim', 'Eq1', 'Ep1', 'EV', 'EP');

%% 测试
swEPSfigure;
figure;
yyaxis left;
plot(beta_nondim, EV, '-ob');
xlabel('$\overline{\beta}$'); ylabel('$E[\overline{V}^2]$');
yyaxis right;
plot(beta_nondim, EP, '--sr');
xlabel('$\overline{\beta}$'); ylabel('$E[\overline{P}]$');
set(legend('$E[\overline{V}^2]$', ...
    '$E[\overline{P}]$', ...
    'Location', 'north'));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22, 'LineWidth', 3.5, 'FontWeight', 'bold');

