%������������ɿ���صĸ�˹���������У����ExΪnumtstep*numbers_noise�ľ���
%tstep�����е�ʱ�䲽��,numtstep�����г���,����������numbers_noise�����������ǿ��Dnȷ��
%������DnΪһά��������һά������ʱ,Ex��������ǿ��������2*(Dn(1),Dn(2),...)������֮������ǿ��Ϊ0
%������DnΪ����ʱ,�Խ�Ԫ��ֵ*2����diag(Dn)*2��Ex�и�������ǿ�ȣ��ǶԽ�Ԫ��*2������֮������ǿ��
%����ͳ�ƽ��cov(ex)��2*Dn/tstep�Ա�
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
Type = 'WhiteNoise'; % Ĭ��ʹ�ø�˹������
[r, c] = size(Dn);
if nargin == 4 && isempty(varargin{4}) == 0
    Type = varargin{4};
end
rr = 0.005; % ����ֲڶ�
v10 = 20; % �ο�����
fs = 10000; % ������
if (r==1 || c==1)
    numbers_noise=length(Dn);
    if norm(Dn)==0; Ex=zeros(numtstep,numbers_noise); return; end  %����Dn��������ʱ
    Ex=zeros(numtstep,numbers_noise);  %����з���������ֵ���з����ǲ�ͬ����
    for k=1:numbers_noise
        s1 = 'WhiteNoise'; s2 = 'Davenport';
        if strcmp(Type, s1) == 1
            Ex(:, k) = randn(numtstep, 1) * sqrt(2 * Dn(k) / tstep); % ��˹������
        elseif strcmp(Type, s2) == 1
            % ����PSD������������
            freqs = (0:numtstep - 1)' / numtstep * fs; % Ƶ������
            f_bar = 1200 * freqs / v10;
            PSD = 4 * rr * v10 ^ 2 * f_bar .^ 2 ./ (freqs .* sign(1 + f_bar .^ 2) .* abs(1 + f_bar .^ 2) .^ (4/3) + eps); % Ŀ�깦�����ܶȣ�����eps�Ա����ĸΪ�㣩
            xx = randn(numtstep, 1); % ��˹������
            w = hamming(numtstep); % ������
            xw = xx .* w; % �Ӵ�
            X = fft(xw); % FFT
            Y = sqrt(PSD) .* X; % Ŀ�깦�����ܶ�
            U_W = real(ifft(Y)); % IFFT
            Ex(:, k) = (U_W / std(U_W)) * sqrt(2 * Dn(k) / tstep); % ��������ǿ��
        else
            error('Error No. 2 in whitenoise.m.');
        end

    end
    return;
end
% if isequal(Dn,Dn')
%     numbers_noise=r;
%     if norm(Dn)==0; Ex=zeros(numtstep,numbers_noise); return; end  %����Dn�������ʱ
%     change=sqrtm(Dn)*sqrt(2/tstep);
%     Ex=randn(numtstep,numbers_noise)*change;  %����з���������ֵ���з����ǲ�ͬ����
%     return;
% end
error('Wrong input parameters in whitenoise.m!');
end