%���������getA()ʹ�ã��õ������ܶ�p�����󣬵���showpxpxy(p,N,domain,C)�ɳ�ͼ
%����4ά����,��Ҫ����p(x2),p(x3),p(x1)��p(x1,x2),p(x1,x3),��
% N=[7,9,11,8];
% Dim=length(N);
% Points=prod(N); domain=rand(Dim,2)+[0,2];
% p=randn(Points,1); C=rand(1,Points);  %����ȡֵ����������ʵ����
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
    if Dim > 1; shij = nchoosek(1:Dim, 2); end %��Ĭ�ϵ�shi��shij��ֵ��nchoosek()����Ϻ���
    if nargin > 6; error('Too many input parameters in showpxpxy !'); end
    if nargin >= 5 && isempty(varargin{5}) == 0; shi = varargin{5}; end %����shi�Ǽ���һά���ݣ���������shi=[2,3,2]����ô���px{1},px{2},px{3}�ֱ��ǵ�2��3��2ά������
    if nargin == 6 && isempty(varargin{6}) == 0; shij = varargin{6}; end %����shij�Ǽ���2ά���ݣ���������shi=[2,3;1,2;1,3]����ô���pxy{1},pxy{2},pxy{3}�ֱ��ǵ�(2,3),(1,2),(1,3)��2ά����
    dx = diff(domain') ./ (N - 1); dV = prod(dx);
    x = cell(1, Dim); for d = 1:Dim; x{d} = linspace(domain(d, 1), domain(d, 2), N(d))'; end %x�Ǹ�ά����
    p = p(:) / (C * p(:) * dV); %���ʹ�һ��
    if Dim == 1; pp = C' .* p * dV; else; pp = reshape(C' .* p * dV, N); end %�Ӹ����ܶȵõ�ÿ����Ԫ�ĸ���
    land = cell(1, Dim);
    for d = 1:Dim; land{d} = ones(N(d), 1); land{d}(1) = 0.5; land{d}(end) = 0.5; end
    toti = length(shi);
    if Dim > 1; totij = size(shij, 1); else; totij = 0; end %����toti��һά���ݣ�totij��2ά����Ҫ����
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

    w = ceil(sqrt(toti + totij)); h = ceil((toti + totij) / w); %��ͼ����w��h��

    if nargout == 0 %���ú���ʱ,���û��������ݣ����ͼ
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
