%本程序产生若干可相关的高斯白噪声序列，输出Ex为numtstep*numbers_noise的矩阵
%tstep是序列的时间步长,numtstep是序列长度,噪声序列数numbers_noise由输入参数半强度Dn确定
%当输入Dn为一维行向量或一维列向量时,Ex各噪声的强度依次是2*(Dn(1),Dn(2),...)，噪声之间的相关强度为0
%当输入Dn为方阵时,对角元素值*2，即diag(Dn)*2是Ex中各噪声的强度，非对角元素*2是噪声之间的相关强度
%例，统计结果cov(ex)与2*Dn/tstep对比
%Dn=[2.5,1.2,0.5;1.2,3.5,1.8;0.5,1.8,2.2];
%tstep=0.02;
%ex=whitenoise(tstep,20000,Dn);
%cov(ex)
%2*Dn/tstep
function Ex = whitenoise(varargin)

if (nargin < 3 || nargin > 4)
    error('Error No. 1 in whitenoise.m.');
end
tstep = varargin{1};
numtstep = varargin{2};
Dn = varargin{3};
Type = 'WhiteNoise'; % 默认使用高斯白噪声
[r, c] = size(Dn);
if nargin == 4 && isempty(varargin{4}) == 0
    Type = varargin{4};
end
rr = 0.005; % 地面粗糙度
v10 = 20; % 参考风速
fs = 10000; % 采样率
if (r==1 || c==1)
    numbers_noise=length(Dn);
    if norm(Dn)==0; Ex=zeros(numtstep,numbers_noise); return; end  %输入Dn是零向量时
    Ex=zeros(numtstep,numbers_noise);  %输出行方向是噪声值，列方向是不同噪声
    for k=1:numbers_noise
        s1 = 'WhiteNoise'; s2 = 'Davenport';
        if strcmp(Type, s1) == 1
            Ex(:, k) = randn(numtstep, 1) * sqrt(2 * Dn(k) / tstep); % 高斯白噪声
        elseif strcmp(Type, s2) == 1
            % 根据PSD生成噪声序列
            freqs = (0:numtstep - 1)' / numtstep * fs; % 频率向量
            f_bar = 1200 * freqs / v10;
            PSD = 4 * rr * v10 ^ 2 * f_bar .^ 2 ./ (freqs .* sign(1 + f_bar .^ 2) .* abs(1 + f_bar .^ 2) .^ (4/3) + eps); % 目标功率谱密度（加了eps以避免分母为零）
            xx = randn(numtstep, 1); % 高斯白噪声
            w = hamming(numtstep); % 汉明窗
            xw = xx .* w; % 加窗
            X = fft(xw); % FFT
            Y = sqrt(PSD) .* X; % 目标功率谱密度
            U_W = real(ifft(Y)); % IFFT
            Ex(:, k) = (U_W / std(U_W)) * sqrt(2 * Dn(k) / tstep); % 调整噪声强度
        else
            error('Error No. 2 in whitenoise.m.');
        end

    end
    return;
end
% if isequal(Dn,Dn')
%     numbers_noise=r;
%     if norm(Dn)==0; Ex=zeros(numtstep,numbers_noise); return; end  %输入Dn是零矩阵时
%     change=sqrtm(Dn)*sqrt(2/tstep);
%     Ex=randn(numtstep,numbers_noise)*change;  %输出行方向是噪声值，列方向是不同噪声
%     return;
% end
error('Wrong input parameters in whitenoise.m!');
end