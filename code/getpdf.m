%输入一维随机向量，统计得到其概率密度函数
%xr,probability_density分别是x轴和概率密度输出，
%例：
%s=randn(1000,1); getpdf(s);
%[x1,p1]=getpdf(s); figure; plot(x1,p1);
%要求第1个参数是一维随机向量
%第2个参数是绘图点数，默认为100，可修改，例如：s=randn(1000,1); getpdf(s,50);
%第3个参数是带宽，默认为自动;
%无输出的调用函数，则绘图，例：s=randn(1000,1); getpdf(s)
function [xr,probability_density] = getpdf(varargin)
samp = varargin{1};
if size(samp,1)~=1 && size(samp,2)~=1; error('Input data of getpdf have to be 1-dimension!' ); end
np = 100;   %默认统计后绘图的点数
Bandwidth = [];    % 默认自动带宽
mini = min(samp); maxi = max(samp); %默认统计区域由样本最小最大值确定
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
[pdf, x] = ksdensity(samp, 'npoints', np, 'Bandwidth', Bandwidth); % 计算平稳概率密度
real_pdf = pdf / sum(pdf(1:end-1) .* diff(x));   % 真实概率密度（归一化）

if nargout==0    %如果无输出，则绘图
    figure;
    plot(x, real_pdf, '.-', 'MarkerSize', 20); xlabel('Element'); ylabel('Probability Density p(Element)');
    return;
end
if nargout==1  %如果输出参数不足，则无数据输出
    return;
end
xr = reshape(x, 1, length(x)); probability_density = reshape(real_pdf, 1, length(real_pdf));
end
