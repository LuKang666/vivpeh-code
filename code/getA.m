function [A, C, ind] = getA(varargin)
    global getA_domain getA_Dim;
    if nargin < 4 || nargin > 5; error('Error No. 1 in getA.'); end %�������Ӧ����4����5��
    class = varargin{end}; N = varargin{end - 1}; getA_domain = varargin{end - 2};
    if ~isvector(N); error('Error No. 2 in getA.'); end %N����������
    N = N(:)'; getA_Dim = length(N); %�����������ΪN-1, ά��getA_Dim
    if ~(ismember(2, size(getA_domain)) && ismember(2, size(class))); error('Error No. 3 in getA.'); end %getA_domain��class��������Ϊ2������Ϊ2
    if size(getA_domain, 2) ~= 2; getA_domain = getA_domain'; end
    if size(class, 2) ~= 2; class = class'; end
    if ~(size(getA_domain, 1) == getA_Dim && size(class, 1) == getA_Dim); error('Error No. 4 in getA.'); end %getA_domain,classά�����붼=getA_Dim
    dx = diff(getA_domain') ./ (N - 1); dV = prod(dx); %��Ԫ����߳�dx�����dV
    Points = prod(N); %������ڵ����������߽��ϵĵ�
    if ismember(0, ismember(class, ["zero", "free", "period"])); error('Error No. 5 in getA.'); end %class��Ԫ��ֻ����"zero","free", "period"�еĴ�
    if ismember(1, sum(ismember(class, "period"), 2)); error('Error No. 6 in getA.'); end %ĳά��"period"�߽�ʱ�������±߽����ͬΪ"period"

    a = zeros(Points, getA_Dim); %����Ư��ϵ������ֵ
    b = zeros(getA_Dim, getA_Dim, Points); %������ɢϵ������ֵ����ɢϵ������Ϊ3*3�����в���Ϊ0��
    getijk = '[ijk(:,1)]=ind2sub(N,1:Points);'; for k = 2:getA_Dim; getijk = insertBefore(getijk, ']=', [',ijk(:,', num2str(k), ')']); end; eval(getijk); %���е��һά���(1:Points)���������ijk
    x = (ijk - 1) .* dx + getA_domain(:, 1)'; %�����е���������������õ�����ֵ��x��ijk����Points*getA_Dim����

    if nargin == 5 %Ư��ϵ��drift����ɢϵ��diffu��getA()�ĵ�һ�Ͷ�����������
        drift = varargin{1}; diffu = varargin{2};
        if length(drift) == 1 && ~iscell(drift); drift = {drift}; end
        if length(diffu) == 1 && ~iscell(diffu); diffu = {diffu}; end
        if ~(iscell(drift) && iscell(diffu)); error('Error No. 7 in getA.'); end %drift��diffu������Ԫ��
        if ~isvector(drift); error('Error No. 8 in getA.'); end %drift����������
        if size(diffu, 1) ~= size(diffu, 2); error('Error No. 9 in getA.'); end %diffu�����Ƿ���
        if ~(length(drift) == 1 || length(drift) == getA_Dim); error('Error No. 10 in getA.'); end %drift����ά��=getA_Dim��=1
        if ~(size(diffu, 1) == 1 || size(diffu, 1) == getA_Dim); error('Error No. 11 in getA.'); end %diffu����ά��=getA_Dim��=1
        if getA_Dim > 1 && length(diffu) ~= 1; for d = 2:getA_Dim; for s = 1:(d - 1); if ~strcmp(diffu{d, s}, "0"); error('Error No. 12 in getA.'); end; end; end; end %diffu�����������Ǿ���������Ԫ��Ϊ"0"���������Խ�Ԫ�أ�

        for no = 1:1 %���������Ư��ϵ��a

            if length(drift) == 1 %���drift�н�һ����������>=2ά��������˵��drift���������ά��Ư�ƺ���ֵ
                if strcmp(drift{1}, "0"); break; end
                type = funtype(drift{1}); %��ȡ��������
                if type == 3; a = feval(drift{1}, x); break; end
                if type == 7; for n = 1:Points; a(n, :) = feval(drift{1}, x(n, :)); end; break; end
            end

            for d = 1:getA_Dim %Ư��ϵ��������
                if strcmp(drift{d}, "0"); continue; end
                type = funtype(drift{d}); %��ȡ��������
                if type == 4; a(:, d) = feval(drift{d}, x); continue; end
                if type == 8 || (type == 6 && getA_Dim == 1) || (type == 7 && getA_Dim == 1); for n = 1:Points; a(n, d) = feval(drift{d}, x(n, :)); end; continue; end
                error('There is function is not being processed.');
            end

        end

        for no = 1:1 %�����������ɢϵ��b

            if size(diffu, 1) == 1 %���diffu��һ����������>=2ά�ķ���˵��diffu���������ά����ɢ����ֵ������������
                if strcmp(diffu{1}, "0"); break; end
                type = funtype(diffu{1});
                if type == 3; b = reshape(feval(diffu{1}, x), 1, 1, Points); break; end
                if type == 2; b = feval(diffu{1}, x); break; end
                if type == 6; for n = 1:Points; b(:, :, n) = feval(diffu{1}, x(n, :)); end; break; end
            end

            for d = 1:getA_Dim

                for s = d:getA_Dim %��ɢϵ���Ǿ���
                    if strcmp(diffu{d, s}, "0"); continue; end
                    type = funtype(diffu{d, s});
                    if type == 4; b(d, s, :) = feval(diffu{d, s}, x); continue; end
                    if type == 8 || (type == 6 && getA_Dim == 1) || (type == 7 && getA_Dim == 1); for n = 1:Points; b(d, s, n) = feval(diffu{d, s}, x(n, :)); end; continue; end
                    error('There is function is not being processed.');
                end

            end

        end

    end

    if nargin == 4 %Ư��ϵ��drift����ɢϵ��diffu����һ����������getA()�ĵ�һ���������
        driftdiffu = varargin{1};
        if length(driftdiffu) ~= 1; error('Error No. 13 in getA.'); end %Ư��ϵ������ɢϵ����1�������м���
        if ~iscell(driftdiffu); driftdiffu = {driftdiffu}; end
        type = funtype(driftdiffu{1}); %��ȡ��������
        if type == 1; [a, b] = feval(driftdiffu{1}, x); end
        if type == 5; for n = 1:Points; [a(n, :), b(:, :, n)] = feval(driftdiffu{1}, x(n, :)); end; end
    end

    drift_yesno = any(a, 1); diffu_yesno = any(b, 3); %��ЩƯ�ƺ���ɢ����ֵΪ�㣨��ȱ��������¼������ȱΪ�㣬��ȱΪ1
    totfun = sum(drift_yesno) + sum(diffu_yesno, [1, 2]); %��¼�¹��ж��ٸ���Ϊ���Ư�ƺ���ɢ���������Ժ�����������ʾ����
    a = a ./ dx; b = b ./ (dx' * dx) ./ (ones(getA_Dim) + eye(getA_Dim)); %ע��Խ�Ԫ��/2�ˣ���Ϊ����ɢ��ǰ��1/2�������

    boundary_1order = {0:1, [-1, 1]; 0:2, [-3, 4, -1] / 2; 0:6, [-147, 360, -450, 400, -225, 72, -10] / 60}; %�߽���һ�׵���ʽ������ֻ�ṩ��3��(�ֱ���һ�׸�ʽ�����׸�ʽ�����׸�ʽ)
    boundary_2order = {0:3, [2, -5, 4, -1]; 0:5, [45, -154, 214, -156, 61, -10] / 12; 0:7, [938, -4014, 7911, -9490, 7380, -3618, 1019, -126] / 180}; %�߽���2�׵���ʽ������ֻ�ṩ��3��(���ף���׺��߽׾��ȸ�ʽ)
    internal_1order = {[-1, 0, 1], [-1, 0, 1] / 2}; %�ڲ�һ�׵���ʽ������ֻ�ṩ��1��
    internal_2order = {[-1, 0, 1], [1, -2, 1]}; %�ڲ�2�׵���ʽ������ֻ�ṩ��1��
    %���ϸ�ʽ��ÿ�ָ�ʽ�������������������ڲ����׵�internal_2order�еĵ�1�֣�internal_2order{1,1}������[-1,0,1]����ʾ���ƽ�����i�ŵ�ʱ��������׵���Ҫ�ĵ��Ϊi+[-1,0,1]��3����ĸ����ܶ�ֵ
    %internal_2order{1,2}������[1,-2,1],������3����Ӧ��Ȩֵ�����Ͽ���⣬��p(j)��ʾһϵ�е��j���ĸ����ܶ�ֵ��f(j)��ʾϵ������ֵ����ôi��f*p��x�Ķ��׵�Ϊp(i+[-1,0,1]).*[1,-2,1].*f(i+[-1,0,1])/dx/dx,���Ƶ�i��f*p��x��һ�׵�Ϊp(i+[-1,0,1]).*[-1,0,1].*f(i+[-1,0,1])/2/dx
    Precision1 = 2; %�������ñ߽���һ�׵������ȣ�������1��3֮�������
    Precision2 = 3; %�������ñ߽��϶��׵������ȣ�������1��3֮�������

    zero_yesno = zeros(1, Points); C = ones(1, Points); % C�ǹ�һ������ p(q1,q2,p1,p2)=Cp(H)

    for d = 1:getA_Dim
        n1 = find(ijk(:, d) == 1); C(n1) = C(n1) / 2; %��¼ÿ������ֵռ��dV�ķ���������p=p/(C*p*dV); �ɶ�p��һ��,pΪ������; n1,n2�ֱ����ÿ��ά���ϵ�����յ�
        if strcmp(class{d, 1}, "zero"); zero_yesno(n1) = 1; end %��¼�õ��Ƿ�����ֵ�߽磬�ǵĻ�Ϊ1,����Ϊ0
        n2 = find(ijk(:, d) == N(d)); C(n2) = C(n2) / 2;
        if strcmp(class{d, 2}, "zero"); zero_yesno(n2) = 1; end
    end

    ind = find(zero_yesno == 0); %�ҵ�������ֵ�߽�ĵ��

    vad = zeros(1, getA_Dim); %��������va,vb��������������ڲ���ʽ���ϱ߽��ʽ���±߽��ʽ������Ҫ����������

    for d = 1:getA_Dim
        if strcmp(class{d, 1}, "period"); vad(d) = 3; continue; end %���ڱ߽�ʱ��va,vb������Ϊ3�м���
        vad(d) = max([3, length(boundary_1order{Precision1, 1}), length(boundary_2order{Precision2, 1})]); %�����ڱ߽�ʱ��va,vb��������ȡΪ�ڲ���ʽ���ϱ߽��ʽ���±߽��ʽ������Ҫ��������
    end

    if ~isempty(find((N - vad) < 0, 1)); error('Error No. 14 in getA.'); end %N��ĳԪ��̫С��Ҫ��Ҫ����va,vb������
    va = cell(1, getA_Dim); vb = cell(1, getA_Dim); %�����ά�����ָ�ʽ��va��һ�׵��ģ�vb��2�׵���
    %va{d}�����ǵ�dάһ�׵��ĸ�ʽ������ΪN(d)����1�����±߽�ĸ�ʽ��Ϊboundary_1order��ĳ1�֣�����N(d)�����ϱ߽�ĸ�ʽ����Ϊboundary_1order��ĳ1�֣�Ҫ�Ӹ��ţ������´���va{d}(N(d),:)=-[...����
    %��2����N(d)-1�����ڲ�����ĸ�ʽ��Ϊinternal_1order�е�ĳ1�֣����߽����ڲ�����Ҫ�ĵ�����һ�£�va{d}���������ȡ���ǵ����ֵ���������߲��㣨�����´���ȵȣ�....[patch,-1,0,1]...��
    for d = 1:getA_Dim

        if strcmp(class(d, 1), "period") %���ڱ߽�
            va{d} = {ones(N(d), 1) * internal_1order{1}, ones(N(d), 1) * internal_1order{2}}; %�������������ڲ��ͱ߽�ĸ�ʽ���߽��ʽһ�ᱻ���ǣ�
            vb{d} = {ones(N(d), 1) * internal_2order{1}, ones(N(d), 1) * internal_2order{2}};
            va{d}{1}(1, :) = [N(d) - 2, 0, 1]; va{d}{2}(1, :) = [-1, 0, 1] / 2; %������4���Ǳ�000���ʽ
            va{d}{1}(N(d), :) = [-N(d) + 2, 0, -1]; va{d}{2}(N(d), :) = [1, 0, -1] / 2;
            vb{d}{1}(1, :) = [N(d) - 2, 0, 1]; vb{d}{2}(1, :) = [1, -2, 1];
            vb{d}{1}(N(d), :) = [-N(d) + 2, 0, -1]; vb{d}{2}(N(d), :) = [1, -2, 1];
            continue;
        end

        if strcmp(class(d, 1), "zero") || strcmp(class(d, 1), "free")
            patch = zeros(1, vad(d) - 3); % ��߲���
            va{d} = {ones(N(d), 1) * [patch, -1, 0, 1], ones(N(d), 1) * [patch, -1, 0, 1] / 2}; %�������������ڲ��ͱ߽�ĸ�ʽ���߽��ʽһ�ᱻ���ǣ�
            vb{d} = {ones(N(d), 1) * [patch, -1, 0, 1], ones(N(d), 1) * [patch, 1, -2, 1]};
            patch = zeros(1, vad(d) - length(boundary_1order{Precision1, 1}));
            va{d}{1}(1, :) = [patch, boundary_1order{Precision1, 1}]; va{d}{2}(1, :) = [patch, boundary_1order{Precision1, 2}]; %������5���Ǳ߽��ʽ
            va{d}{1}(N(d), :) = -va{d}{1}(1, :); va{d}{2}(N(d), :) = -va{d}{2}(1, :);
            patch = zeros(1, vad(d) - length(boundary_2order{Precision2, 1}));
            vb{d}{1}(1, :) = [patch, boundary_2order{Precision2, 1}]; vb{d}{2}(1, :) = [patch, boundary_2order{Precision2, 2}];
            vb{d}{1}(N(d), :) = -vb{d}{1}(1, :); vb{d}{2}(N(d), :) = vb{d}{2}(1, :);
        end

    end

    e = eye(getA_Dim);
    totf = 0;
    if Points > 500; A = sparse(Points, Points); Ax = sparse(Points, Points); else A = zeros(Points); Ax = zeros(Points); end %����500�����Ӧ��ϡ������

    if getA_Dim == 1
        run1 = 'reshape(no(:,1,:),Points,vad(d));';
    else
        run1 = 'sub2ind(N);'; for k = 1:getA_Dim; run1 = insertBefore(run1, ');', [',reshape(no(:,', num2str(k), ',:),Points,vad(d))']); end
        run2 = 'sub2ind(N);'; for k = 1:getA_Dim; run2 = insertBefore(run2, ');', [',reshape(no(:,', num2str(k), ',:),Points,vad(d)*vad(s))']); end
    end

    for d = 1:getA_Dim

        if drift_yesno(d) %Ư��������ڣ�������
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            vv = va{d}{1}(ijk(:, d), :);
            no = ijk + permute(reshape(vv(:) * e(d, :), Points, vad(d), getA_Dim), [1, 3, 2]);
            no2 = eval(run1);
            cc = -va{d}{2}(ijk(:, d), :);
            z = (1:Points)' * ones(1, vad(d));
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            Ax = spconvert([z(:), no2(:), cc(:); Points, Points, 0]);
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            A = A + Ax .* a(:, d).'; % ���¸����ܶȾ���
        end

        if diffu_yesno(d, d) %��ɢ������ڣ�������
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            vv = vb{d}{1}(ijk(:, d), :);
            no = ijk + permute(reshape(vv(:) * e(d, :), Points, vad(d), getA_Dim), [1, 3, 2]);
            no2 = eval(run1);
            cc = vb{d}{2}(ijk(:, d), :);
            z = (1:Points)' * ones(1, vad(d));
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            Ax = spconvert([z(:), no2(:), cc(:); Points, Points, 0]);
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            A = A + Ax .* reshape(b(d, d, :), [1, Points]);
        end

        for s = (d + 1):getA_Dim

            if diffu_yesno(d, s) %����������ڣ�������
                totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
                v1 = va{d}{1}(ijk(:, d), :);
                v2 = va{s}{1}(ijk(:, s), :)';
                v1 = v1(:) * ones(1, vad(s));
                v2 = reshape(v2(:) * ones(1, vad(d)), vad(s), Points * vad(d))';
                no = ijk + permute(reshape(v1(:) * e(d, :) + v2(:) * e(s, :), Points, vad(d) * vad(s), getA_Dim), [1, 3, 2]);
                totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
                no2 = eval(run2);
                w1 = va{d}{2}(ijk(:, d), :);
                w1 = w1(:) * ones(1, vad(s));
                w2 = va{s}{2}(ijk(:, s), :)';
                w2 = reshape(w2(:) * ones(1, vad(d)), vad(s), Points * vad(d))';
                cc = reshape(w1(:) .* w2(:), Points, vad(d) * vad(s));
                z = (1:Points)' * ones(1, vad(d) * vad(s));
                totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
                Ax = spconvert([z(:), no2(:), cc(:); Points, Points, 0]);
                A = A + Ax .* reshape(b(d, s, :), [1, Points]);
            end

        end

    end

    disp('Finished getA.');
    return;
end

function type = funtype(fun) %���ݺ�����������ж�����
    global getA_domain getA_Dim;
    type = 0;
    x = (0:20)' * diff(getA_domain') / 20 + getA_domain(:, 1)'; %x��[21,getA_Dim]ά�Դ�����
    y = []; y1 = []; y2 = [];

    for no = 1:1
        try [y1, y2] = feval(fun, x); break; end
        try [y1, y2] = feval(fun, x(10, :)); break; end
        try y = feval(fun, x); break; end
        try y = feval(fun, x(10, :)); break; end
    end

    for no = 1:1
        try if ~(isequal(size(y1), [21, getA_Dim]) && isequal(size(y2), [getA_Dim, getA_Dim, 21])); error('*'); end; type = 1; break; end %1��ָͬʱ����[Points,getA_Dim]ά��[getA_Dim,getA_Dim,Points]ά���ݵĺ�������ͬʱ��������Ư��ϵ����������ɢϵ���ĺ�����
        try if ~isequal(size(y), [getA_Dim, getA_Dim, 21]); error('*'); end; type = 2; break; end %2��ָ����[getA_Dim,getA_Dim,Points]ά���ݵĺ������ܼ���������ɢϵ���ĺ�����
        try if ~isequal(size(y), [21, getA_Dim]); error('*'); end; type = 3; break; end %3��ָ����[Points,getA_Dim]ά���ݵĺ������ܼ�������Ư��ϵ���ĺ�����
        try if ~(isequal(size(y), [21, 1]) | isequal(size(y), [1, 21])); error('*'); end; type = 4; break; end %4��ָ����[Points,1]ά����[1,Points]ά�����ݵĺ������ܼ���Ư��ϵ���������ɢϵ���ĺ�����
        try if ~((isequal(size(y1), [1, getA_Dim]) | isequal(size(y1), [getA_Dim, 1])) && isequal(size(y2), [getA_Dim, getA_Dim])); error('*'); end; type = 5; break; end %5��ָͬʱ����[getA_Dim,1]ά����[1,getA_Dim]ά����[getA_Dim,getA_Dim]ά���ݵĺ�������ͬʱ��������Ư��ϵ����������ɢϵ���ĺ�����
        try if ~isequal(size(y), [getA_Dim, getA_Dim]); error('*'); end; type = 6; break; end %6��ָ����[getA_Dim,getA_Dim]ά���ݵĺ������ܼ���������ɢϵ���ĺ�����
        try if ~(isequal(size(y), [1, getA_Dim]) | isequal(size(y), [getA_Dim, 1])); error('*'); end; type = 7; break; end %7��ָ����[getA_Dim,1]ά����[1,getA_Dim]ά�����ݵĺ������ܼ�������Ư��ϵ���ĺ�����
        try if ~isequal(size(y), [1, 1]); error('*'); end; type = 8; break; end %8��ָ����[1,1]ά���ݵĺ������ܼ���1��Ư��ϵ����1����ɢϵ���ĺ�����
    end

    if type == 0; error('The input function has no corresponding type.'); end
end
