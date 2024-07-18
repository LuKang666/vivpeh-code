%本程序配合getA()使用，得到概率密度p向量后，调用showpxpxy(p,N,domain,C)可出图
%比如4维情形,若要导出p(x2),p(x3),p(x1)和p(x1,x2),p(x1,x3),则：
% N=[7,9,11,8];
% Dim=length(N);
% Points=prod(N); domain=rand(Dim,2)+[0,2];
% p=randn(Points,1); C=rand(1,Points);  %任意取值，不代表真实数据
% showpxpxy(p,N,domain,C);
% [x,px,pxy]=showpxpxy(p,N,domain,C,[2,3,1],[1,2;1,3]);
% figure; plot(x{2},px{1});
% figure; plot(x{3},px{2});
% figure; plot(x{1},px{3});
% figure; surf(x{1},x{2},pxy{1}');
% figure; surf(x{1},x{3},pxy{2}');
function [x, px, pxy] = showpxpxy(varargin)
    p = varargin{1}; N = varargin{2}; domain = varargin{3}; C = varargin{4};
    Dim = length(N); shi = 1:Dim;
    if Dim > 1; shij = nchoosek(1:Dim, 2); end %对默认的shi和shij赋值，nchoosek()是组合函数
    if nargin > 6; error('Too many input parameters in showpxpxy !'); end
    if nargin >= 5 && isempty(varargin{5}) == 0; shi = varargin{5}; end %输入shi是计算一维数据，比如输入shi=[2,3,2]，那么输出px{1},px{2},px{3}分别是第2，3，2维的数据
    if nargin == 6 && isempty(varargin{6}) == 0; shij = varargin{6}; end %输入shij是计算2维数据，比如输入shi=[2,3;1,2;1,3]，那么输出pxy{1},pxy{2},pxy{3}分别是第(2,3),(1,2),(1,3)的2维数据
    dx = diff(domain') ./ (N - 1); dV = prod(dx);
    x = cell(1, Dim); for d = 1:Dim; x{d} = linspace(domain(d, 1), domain(d, 2), N(d))'; end %x是各维坐标
    p = p(:) / (C * p(:) * dV); %概率归一化
    if Dim == 1; pp = C' .* p * dV; else; pp = reshape(C' .* p * dV, N); end %从概率密度得到每个单元的概率
    land = cell(1, Dim);
    for d = 1:Dim; land{d} = ones(N(d), 1); land{d}(1) = 0.5; land{d}(end) = 0.5; end
    toti = length(shi);
    if Dim > 1; totij = size(shij, 1); else; totij = 0; end %共有toti个一维数据，totij个2维数据要计算
    px = cell(1, toti);
    if Dim > 1; pxy = cell(1, totij); end
    if Dim == 1; pxy = NaN; end

    for k = 1:toti
        i = shi(k);
        if Dim == 1; px{k} = pp ./ land{i} / dx(i); else; px{k} = reshape(sum(pp, setdiff(1:Dim, i)), [N(i), 1]) ./ land{i} / dx(i); end
    end

    if Dim > 1

        for k = 1:totij
            i = shij(k, 1); j = shij(k, 2);
            integration = setdiff(1:Dim, [i, j]);

            if isempty(integration)
                pxy{k} = reshape(pp, N([i, j])) ./ (land{i} * land{j}') / dx(i) / dx(j);
            else
                pxy{k} = reshape(sum(pp, integration), N([i, j])) ./ (land{i} * land{j}') / dx(i) / dx(j);
            end

        end

    end

    w = ceil(sqrt(toti + totij)); h = ceil((toti + totij) / w); %绘图安排w列h行

    if nargout == 0 %调用函数时,如果没有输出数据，则绘图
        figure;

        for k = 1:toti
            i = shi(k); subplot(h, w, k); plot(x{i}, px{k}, '-k'); xlabel(['H', num2str(i)]);
        end

        if Dim > 1

            for k = 1:totij
                i = shij(k, 1); j = shij(k, 2); subplot(h, w, k + toti); surf(x{i}, x{j}, pxy{k}'); xlabel(['H', num2str(i)]); ylabel(['H', num2str(j)]);
            end

        end

    end

end
