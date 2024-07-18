%输入一维随机向量，统计得到其概率密度函数
%xr,probability_density分别是x轴和概率密度输出，
%例：
%s=randn(1000,1); getpdf(s);
%[x1,p1]=getpdf(s); figure; plot(x1,p1);
%要求第1个参数是一维随机向量
%第2，3个参数是统计的下，上限，例：s=randn(1000,1); getpdf(s,-2,1.5);
%缺省第2，3参数或输入空值[]时，则从随机向量中自动计算最小最大值
%第4个参数是绘图点数，默认为20，可修改，例如：s=randn(1000,1); getpdf(s,[],[],50);
%无输出的调用函数，则绘图，例：s=randn(1000,1); getpdf(s)
function [xr, probability_density] = getpdf_1(varargin)
    samp = varargin{1};
    if size(samp, 1) ~= 1 && size(samp, 2) ~= 1; error('Input data of getpdf have to be 1-dimension!'); end
    np = 20; %默认统计后绘图的点数
    mini = min(samp); maxi = max(samp); %默认统计区域由样本最小最大值确定
    if nargin > 3; error('Too many input parameters in getpdf!'); end

    if nargin >= 2 && isempty(varargin{2}) == 0
        range = varargin{2};
        mini = range(1);
        maxi = range(2);
    end

    if nargin >= 3 && isempty(varargin{3}) == 0
        np = varargin{3};
    end


    Xedges = linspace(mini, maxi, np + 1); %确定间隔
    x = (Xedges(1:end - 1) + Xedges(2:end)) / 2; %取间隔的中间值做为x轴的输出
    p = histcounts(samp, Xedges) / length(samp) / (Xedges(2) - Xedges(1)); %histcounts函数是统计落入间隔内的数目

    if nargout == 0 %如果无输出，则绘图
        figure;
        plot(x, p, '.-', 'MarkerSize', 20); xlabel('Element'); ylabel('Probability Density p(Element)');
        return;
    end

    if nargout == 1 %如果输出参数不足，则无数据输出
        return;
    end

    xr = x'; probability_density = p';
end
