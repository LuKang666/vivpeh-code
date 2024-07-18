%����һά���������ͳ�Ƶõ�������ܶȺ���
%xr,probability_density�ֱ���x��͸����ܶ������
%����
%s=randn(1000,1); getpdf(s);
%[x1,p1]=getpdf(s); figure; plot(x1,p1);
%Ҫ���1��������һά�������
%��2�������ǻ�ͼ������Ĭ��Ϊ100�����޸ģ����磺s=randn(1000,1); getpdf(s,50);
%��3�������Ǵ���Ĭ��Ϊ�Զ�;
%������ĵ��ú��������ͼ������s=randn(1000,1); getpdf(s)
function [xr,probability_density] = getpdf(varargin)
samp = varargin{1};
if size(samp,1)~=1 && size(samp,2)~=1; error('Input data of getpdf have to be 1-dimension!' ); end
np = 100;   %Ĭ��ͳ�ƺ��ͼ�ĵ���
Bandwidth = [];    % Ĭ���Զ�����
mini = min(samp); maxi = max(samp); %Ĭ��ͳ��������������С���ֵȷ��
if nargin > 4; error('Too many input parameters in getpdf!'); end
if nargin >= 2 && isempty(varargin{2}) == 0
    range = varargin{2};
    mini = range(1);
    maxi = range(2);
end

if nargin >= 3 && isempty(varargin{3}) == 0
    np = varargin{3};
end

if nargin == 4 && isempty(varargin{4}) == 0
    Bandwidth = varargin{4};
end


samp = linspace(mini, maxi, np + 1);
% [pdf, x] = ksdensity(samp, 'npoints', np);
[pdf, x] = ksdensity(samp, 'npoints', np, 'Bandwidth', Bandwidth); % ����ƽ�ȸ����ܶ�
real_pdf = pdf / sum(pdf(1:end-1) .* diff(x));   % ��ʵ�����ܶȣ���һ����

if nargout==0    %�������������ͼ
    figure;
    plot(x, real_pdf, '.-', 'MarkerSize', 20); xlabel('Element'); ylabel('Probability Density p(Element)');
    return;
end
if nargout==1  %�������������㣬�����������
    return;
end
xr = reshape(x, 1, length(x)); probability_density = reshape(real_pdf, 1, length(real_pdf));
end
