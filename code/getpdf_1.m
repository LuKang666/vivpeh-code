%����һά���������ͳ�Ƶõ�������ܶȺ���
%xr,probability_density�ֱ���x��͸����ܶ������
%����
%s=randn(1000,1); getpdf(s);
%[x1,p1]=getpdf(s); figure; plot(x1,p1);
%Ҫ���1��������һά�������
%��2��3��������ͳ�Ƶ��£����ޣ�����s=randn(1000,1); getpdf(s,-2,1.5);
%ȱʡ��2��3�����������ֵ[]ʱ���������������Զ�������С���ֵ
%��4�������ǻ�ͼ������Ĭ��Ϊ20�����޸ģ����磺s=randn(1000,1); getpdf(s,[],[],50);
%������ĵ��ú��������ͼ������s=randn(1000,1); getpdf(s)
function [xr, probability_density] = getpdf_1(varargin)
    samp = varargin{1};
    if size(samp, 1) ~= 1 && size(samp, 2) ~= 1; error('Input data of getpdf have to be 1-dimension!'); end
    np = 20; %Ĭ��ͳ�ƺ��ͼ�ĵ���
    mini = min(samp); maxi = max(samp); %Ĭ��ͳ��������������С���ֵȷ��
    if nargin > 3; error('Too many input parameters in getpdf!'); end

    if nargin >= 2 && isempty(varargin{2}) == 0
        range = varargin{2};
        mini = range(1);
        maxi = range(2);
    end

    if nargin >= 3 && isempty(varargin{3}) == 0
        np = varargin{3};
    end


    Xedges = linspace(mini, maxi, np + 1); %ȷ�����
    x = (Xedges(1:end - 1) + Xedges(2:end)) / 2; %ȡ������м�ֵ��Ϊx������
    p = histcounts(samp, Xedges) / length(samp) / (Xedges(2) - Xedges(1)); %histcounts������ͳ���������ڵ���Ŀ

    if nargout == 0 %�������������ͼ
        figure;
        plot(x, p, '.-', 'MarkerSize', 20); xlabel('Element'); ylabel('Probability Density p(Element)');
        return;
    end

    if nargout == 1 %�������������㣬�����������
        return;
    end

    xr = x'; probability_density = p';
end
