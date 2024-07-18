%输入2维随机向量，统计得到其概率密度函数
%xr,yr,probability_density分别是x,y轴和概率密度输出，例
%s=[randn(1000,1),rand(1000,1)]; [x,y,p]=getpdf2(s); surf(x,y,p'); xlabel('x'); ylabel('y');
%要求第1个参数是2维随机向量矩阵，行数为2（则第1行是x值，第2行是y值）或列数为2（则第1列是x值，第2列是y值）都可以
%第2，3个参数分别是x,y统计的下，上限
%例：s=[randn(1000,1),rand(1000,1)]; getpdf2(s,[-3,3],[0,1]);表明第一维统计范围-3至3，第2维0至1
%缺省第2，3参数时，则自动取samp的最小最大值
%第4个参数是最小间隔，默认是0.1*0.1，可修改，例改成0.01*0.02：s=[randn(1000,1),rand(1000,1)]; getpdf2(s,[],[],[0.01,0.02])
%第5个参数是带宽，默认为自动;
%无输出的调用函数，则绘图，例：s=[randn(1000,1),rand(1000,1)]; getpdf2(s)
function [xr,yr,probability_density]=getpdf2(varargin)
samp=varargin{1};
if size(samp,1)<size(samp,2); samp=samp'; end
if size(samp,2)~=2; error('Input data of getpdf2 have to be 2-dimension!' ); end
npx = 20; npy = 20; %默认统计后绘图的点数
mini=min(samp); maxi=max(samp);
Bandwidth = []; % 默认自动带宽

if nargin>5; error('Too many input parameters in getpdf2!'); end
if nargin>=2 && isempty(varargin{2})==0
    mini(1)=varargin{2}(1);
    maxi(1)=varargin{2}(2);
end
if nargin>=3 && isempty(varargin{3})==0
    mini(2)=varargin{3}(1);
    maxi(2)=varargin{3}(2);
end
if nargin>=4 && isempty(varargin{4})==0
    npx = varargin{4}(1);
    npy = varargin{4}(2);
end
if nargin == 5 && isempty(varargin{5}) == 0
    Bandwidth = varargin{5};
end

x_range = linspace(mini(1), maxi(1), npx);  % 确定x轴间隔
y_range = linspace(mini(2), maxi(2), npy);  % 确定y轴间隔
[X, Y] = meshgrid(x_range, y_range);
XY = [X(:), Y(:)];
pdf = mvksdensity(samp, XY, 'Bandwidth', Bandwidth); % 计算联合概率密度
dx = (maxi(1) - mini(1)) / (npx - 1);
dy = (maxi(2) - mini(2)) / (npy - 1);
real_pdf = pdf / sum(pdf * dx * dy, 'all'); % 真实概率密度(归一化)

if nargout==0    %如果无输出，则绘图
    figure;
    surf(X,Y,reshape(real_pdf,size(X))); xlabel('Element 1'); ylabel('Element 2'); zlabel('Probability Density p(Element 1, Element 2)');
    return;
end
if nargout==1 || nargout==2  %如果输出参数不足，则无数据输出
    return;
end
xr = X; yr = Y; probability_density = reshape(real_pdf, size(X));
end