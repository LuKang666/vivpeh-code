%����2ά���������ͳ�Ƶõ�������ܶȺ���
%xr,yr,probability_density�ֱ���x,y��͸����ܶ��������
%s=[randn(1000,1),rand(1000,1)]; [x,y,p]=getpdf2(s); surf(x,y,p'); xlabel('x'); ylabel('y');
%Ҫ���1��������2ά���������������Ϊ2�����1����xֵ����2����yֵ��������Ϊ2�����1����xֵ����2����yֵ��������
%��2��3�������ֱ���x,yͳ�Ƶ��£�����
%����s=[randn(1000,1),rand(1000,1)]; getpdf2(s,[-3,3],[0,1]);������һάͳ�Ʒ�Χ-3��3����2ά0��1
%ȱʡ��2��3����ʱ�����Զ�ȡsamp����С���ֵ
%��4����������С�����Ĭ����0.1*0.1�����޸ģ����ĳ�0.01*0.02��s=[randn(1000,1),rand(1000,1)]; getpdf2(s,[],[],[0.01,0.02])
%��5�������Ǵ���Ĭ��Ϊ�Զ�;
%������ĵ��ú��������ͼ������s=[randn(1000,1),rand(1000,1)]; getpdf2(s)
function [xr,yr,probability_density]=getpdf2(varargin)
samp=varargin{1};
if size(samp,1)<size(samp,2); samp=samp'; end
if size(samp,2)~=2; error('Input data of getpdf2 have to be 2-dimension!' ); end
npx = 20; npy = 20; %Ĭ��ͳ�ƺ��ͼ�ĵ���
mini=min(samp); maxi=max(samp);
Bandwidth = []; % Ĭ���Զ�����

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

x_range = linspace(mini(1), maxi(1), npx);  % ȷ��x����
y_range = linspace(mini(2), maxi(2), npy);  % ȷ��y����
[X, Y] = meshgrid(x_range, y_range);
XY = [X(:), Y(:)];
pdf = mvksdensity(samp, XY, 'Bandwidth', Bandwidth); % �������ϸ����ܶ�
dx = (maxi(1) - mini(1)) / (npx - 1);
dy = (maxi(2) - mini(2)) / (npy - 1);
real_pdf = pdf / sum(pdf * dx * dy, 'all'); % ��ʵ�����ܶ�(��һ��)

if nargout==0    %�������������ͼ
    figure;
    surf(X,Y,reshape(real_pdf,size(X))); xlabel('Element 1'); ylabel('Element 2'); zlabel('Probability Density p(Element 1, Element 2)');
    return;
end
if nargout==1 || nargout==2  %�������������㣬�����������
    return;
end
xr = X; yr = Y; probability_density = reshape(real_pdf, size(X));
end