function [A, C, ind] = getA(varargin)
    global getA_domain getA_Dim;
    if nargin < 4 || nargin > 5; error('Error No. 1 in getA.'); end %输入参数应该是4个或5个
    class = varargin{end}; N = varargin{end - 1}; getA_domain = varargin{end - 2};
    if ~isvector(N); error('Error No. 2 in getA.'); end %N必须是向量
    N = N(:)'; getA_Dim = length(N); %求解域网格数为N-1, 维数getA_Dim
    if ~(ismember(2, size(getA_domain)) && ismember(2, size(class))); error('Error No. 3 in getA.'); end %getA_domain和class必须行数为2或列数为2
    if size(getA_domain, 2) ~= 2; getA_domain = getA_domain'; end
    if size(class, 2) ~= 2; class = class'; end
    if ~(size(getA_domain, 1) == getA_Dim && size(class, 1) == getA_Dim); error('Error No. 4 in getA.'); end %getA_domain,class维数必须都=getA_Dim
    dx = diff(getA_domain') ./ (N - 1); dV = prod(dx); %单元体各边长dx和体积dV
    Points = prod(N); %求解域内点数，包括边界上的点
    if ismember(0, ismember(class, ["zero", "free", "period"])); error('Error No. 5 in getA.'); end %class的元素只能是"zero","free", "period"中的词
    if ismember(1, sum(ismember(class, "period"), 2)); error('Error No. 6 in getA.'); end %某维是"period"边界时，其上下边界必须同为"period"

    a = zeros(Points, getA_Dim); %存贮漂移系数函数值
    b = zeros(getA_Dim, getA_Dim, Points); %存贮扩散系数函数值（扩散系数个数为3*3，其中部分为0）
    getijk = '[ijk(:,1)]=ind2sub(N,1:Points);'; for k = 2:getA_Dim; getijk = insertBefore(getijk, ']=', [',ijk(:,', num2str(k), ')']); end; eval(getijk); %所有点的一维编号(1:Points)换成坐标号ijk
    x = (ijk - 1) .* dx + getA_domain(:, 1)'; %从所有点的坐标线性索引得到坐标值（x和ijk都是Points*getA_Dim矩阵）

    if nargin == 5 %漂移系数drift和扩散系数diffu在getA()的第一和二个参数输入
        drift = varargin{1}; diffu = varargin{2};
        if length(drift) == 1 && ~iscell(drift); drift = {drift}; end
        if length(diffu) == 1 && ~iscell(diffu); diffu = {diffu}; end
        if ~(iscell(drift) && iscell(diffu)); error('Error No. 7 in getA.'); end %drift和diffu必须是元胞
        if ~isvector(drift); error('Error No. 8 in getA.'); end %drift必须是向量
        if size(diffu, 1) ~= size(diffu, 2); error('Error No. 9 in getA.'); end %diffu必须是方阵
        if ~(length(drift) == 1 || length(drift) == getA_Dim); error('Error No. 10 in getA.'); end %drift必须维数=getA_Dim或=1
        if ~(size(diffu, 1) == 1 || size(diffu, 1) == getA_Dim); error('Error No. 11 in getA.'); end %diffu必须维数=getA_Dim或=1
        if getA_Dim > 1 && length(diffu) ~= 1; for d = 2:getA_Dim; for s = 1:(d - 1); if ~strcmp(diffu{d, s}, "0"); error('Error No. 12 in getA.'); end; end; end; end %diffu必须是上三角矩阵，下三角元素为"0"（不包括对角元素）

        for no = 1:1 %计算和贮存漂移系数a

            if length(drift) == 1 %如果drift中仅一个函数而非>=2维的向量，说明drift能输出所有维的漂移函数值
                if strcmp(drift{1}, "0"); break; end
                type = funtype(drift{1}); %获取函数类型
                if type == 3; a = feval(drift{1}, x); break; end
                if type == 7; for n = 1:Points; a(n, :) = feval(drift{1}, x(n, :)); end; break; end
            end

            for d = 1:getA_Dim %漂移系数是向量
                if strcmp(drift{d}, "0"); continue; end
                type = funtype(drift{d}); %获取函数类型
                if type == 4; a(:, d) = feval(drift{d}, x); continue; end
                if type == 8 || (type == 6 && getA_Dim == 1) || (type == 7 && getA_Dim == 1); for n = 1:Points; a(n, d) = feval(drift{d}, x(n, :)); end; continue; end
                error('There is function is not being processed.');
            end

        end

        for no = 1:1 %计算和贮存扩散系数b

            if size(diffu, 1) == 1 %如果diffu仅一个函数而非>=2维的方阵，说明diffu能输出所有维的扩散函数值，包括交叉项
                if strcmp(diffu{1}, "0"); break; end
                type = funtype(diffu{1});
                if type == 3; b = reshape(feval(diffu{1}, x), 1, 1, Points); break; end
                if type == 2; b = feval(diffu{1}, x); break; end
                if type == 6; for n = 1:Points; b(:, :, n) = feval(diffu{1}, x(n, :)); end; break; end
            end

            for d = 1:getA_Dim

                for s = d:getA_Dim %扩散系数是矩阵
                    if strcmp(diffu{d, s}, "0"); continue; end
                    type = funtype(diffu{d, s});
                    if type == 4; b(d, s, :) = feval(diffu{d, s}, x); continue; end
                    if type == 8 || (type == 6 && getA_Dim == 1) || (type == 7 && getA_Dim == 1); for n = 1:Points; b(d, s, n) = feval(diffu{d, s}, x(n, :)); end; continue; end
                    error('There is function is not being processed.');
                end

            end

        end

    end

    if nargin == 4 %漂移系数drift和扩散系数diffu共用一个函数，在getA()的第一个输入参数
        driftdiffu = varargin{1};
        if length(driftdiffu) ~= 1; error('Error No. 13 in getA.'); end %漂移系数和扩散系数在1个函数中计算
        if ~iscell(driftdiffu); driftdiffu = {driftdiffu}; end
        type = funtype(driftdiffu{1}); %获取函数类型
        if type == 1; [a, b] = feval(driftdiffu{1}, x); end
        if type == 5; for n = 1:Points; [a(n, :), b(:, :, n)] = feval(driftdiffu{1}, x(n, :)); end; end
    end

    drift_yesno = any(a, 1); diffu_yesno = any(b, 3); %有些漂移和扩散函数值为零（或缺乏），记录下来，缺为零，不缺为1
    totfun = sum(drift_yesno) + sum(diffu_yesno, [1, 2]); %记录下共有多少个不为零的漂移和扩散函数，用以后续计算中显示进度
    a = a ./ dx; b = b ./ (dx' * dx) ./ (ones(getA_Dim) + eye(getA_Dim)); %注意对角元素/2了，因为把扩散项前的1/2搬进来了

    boundary_1order = {0:1, [-1, 1]; 0:2, [-3, 4, -1] / 2; 0:6, [-147, 360, -450, 400, -225, 72, -10] / 60}; %边界上一阶导格式，这里只提供了3种(分别是一阶格式，二阶格式和六阶格式)
    boundary_2order = {0:3, [2, -5, 4, -1]; 0:5, [45, -154, 214, -156, 61, -10] / 12; 0:7, [938, -4014, 7911, -9490, 7380, -3618, 1019, -126] / 180}; %边界上2阶导格式，这里只提供了3种(二阶，五阶和七阶精度格式)
    internal_1order = {[-1, 0, 1], [-1, 0, 1] / 2}; %内部一阶导格式，这里只提供了1种
    internal_2order = {[-1, 0, 1], [1, -2, 1]}; %内部2阶导格式，这里只提供了1种
    %以上格式中每种格式都有两个向量，比如内部二阶导internal_2order中的第1种，internal_2order{1,1}是向量[-1,0,1]，表示当推进至第i号点时，计算二阶导需要的点号为i+[-1,0,1]共3个点的概率密度值
    %internal_2order{1,2}是向量[1,-2,1],是上述3点相应的权值，综上可理解，令p(j)表示一系列点号j处的概率密度值，f(j)表示系数函数值，那么i处f*p对x的二阶导为p(i+[-1,0,1]).*[1,-2,1].*f(i+[-1,0,1])/dx/dx,类似的i处f*p对x的一阶导为p(i+[-1,0,1]).*[-1,0,1].*f(i+[-1,0,1])/2/dx
    Precision1 = 2; %这里设置边界上一阶导数精度，必须是1至3之间的整数
    Precision2 = 3; %这里设置边界上二阶导数精度，必须是1至3之间的整数

    zero_yesno = zeros(1, Points); C = ones(1, Points); % C是归一化常数 p(q1,q2,p1,p2)=Cp(H)

    for d = 1:getA_Dim
        n1 = find(ijk(:, d) == 1); C(n1) = C(n1) / 2; %记录每个概率值占用dV的份数，运行p=p/(C*p*dV); 可对p归一化,p为列向量; n1,n2分别代表每个维度上的起点终点
        if strcmp(class{d, 1}, "zero"); zero_yesno(n1) = 1; end %记录该点是否在零值边界，是的话为1,不是为0
        n2 = find(ijk(:, d) == N(d)); C(n2) = C(n2) / 2;
        if strcmp(class{d, 2}, "zero"); zero_yesno(n2) = 1; end
    end

    ind = find(zero_yesno == 0); %找到不是零值边界的点号

    vad = zeros(1, getA_Dim); %贮存以下va,vb矩阵的列数（即内部格式，上边界格式，下边界格式中所需要的最多点数）

    for d = 1:getA_Dim
        if strcmp(class{d, 1}, "period"); vad(d) = 3; continue; end %周期边界时，va,vb矩阵定义为3列即可
        vad(d) = max([3, length(boundary_1order{Precision1, 1}), length(boundary_2order{Precision2, 1})]); %非周期边界时，va,vb矩阵列数取为内部格式，上边界格式，下边界格式中所需要的最多点数
    end

    if ~isempty(find((N - vad) < 0, 1)); error('Error No. 14 in getA.'); end %N的某元素太小，要求要大于va,vb的列数
    va = cell(1, getA_Dim); vb = cell(1, getA_Dim); %贮存各维网格差分格式，va是一阶导的，vb是2阶导的
    %va{d}矩阵是第d维一阶导的格式，行数为N(d)，第1行是下边界的格式（为boundary_1order的某1种），第N(d)行是上边界的格式（亦为boundary_1order的某1种，要加负号，见以下代码va{d}(N(d),:)=-[...），
    %第2行至N(d)-1行是内部区域的格式（为internal_1order中的某1种），边界与内部所需要的点数不一致，va{d}矩阵的列数取他们的最大值，不足的左边补零（见以下代码等等，....[patch,-1,0,1]...）
    for d = 1:getA_Dim

        if strcmp(class(d, 1), "period") %周期边界
            va{d} = {ones(N(d), 1) * internal_1order{1}, ones(N(d), 1) * internal_1order{2}}; %这以下两行是内部和边界的格式（边界格式一会被覆盖）
            vb{d} = {ones(N(d), 1) * internal_2order{1}, ones(N(d), 1) * internal_2order{2}};
            va{d}{1}(1, :) = [N(d) - 2, 0, 1]; va{d}{2}(1, :) = [-1, 0, 1] / 2; %这以下4行是边000界格式
            va{d}{1}(N(d), :) = [-N(d) + 2, 0, -1]; va{d}{2}(N(d), :) = [1, 0, -1] / 2;
            vb{d}{1}(1, :) = [N(d) - 2, 0, 1]; vb{d}{2}(1, :) = [1, -2, 1];
            vb{d}{1}(N(d), :) = [-N(d) + 2, 0, -1]; vb{d}{2}(N(d), :) = [1, -2, 1];
            continue;
        end

        if strcmp(class(d, 1), "zero") || strcmp(class(d, 1), "free")
            patch = zeros(1, vad(d) - 3); % 左边补零
            va{d} = {ones(N(d), 1) * [patch, -1, 0, 1], ones(N(d), 1) * [patch, -1, 0, 1] / 2}; %这以下两行是内部和边界的格式（边界格式一会被覆盖）
            vb{d} = {ones(N(d), 1) * [patch, -1, 0, 1], ones(N(d), 1) * [patch, 1, -2, 1]};
            patch = zeros(1, vad(d) - length(boundary_1order{Precision1, 1}));
            va{d}{1}(1, :) = [patch, boundary_1order{Precision1, 1}]; va{d}{2}(1, :) = [patch, boundary_1order{Precision1, 2}]; %这以下5行是边界格式
            va{d}{1}(N(d), :) = -va{d}{1}(1, :); va{d}{2}(N(d), :) = -va{d}{2}(1, :);
            patch = zeros(1, vad(d) - length(boundary_2order{Precision2, 1}));
            vb{d}{1}(1, :) = [patch, boundary_2order{Precision2, 1}]; vb{d}{2}(1, :) = [patch, boundary_2order{Precision2, 2}];
            vb{d}{1}(N(d), :) = -vb{d}{1}(1, :); vb{d}{2}(N(d), :) = vb{d}{2}(1, :);
        end

    end

    e = eye(getA_Dim);
    totf = 0;
    if Points > 500; A = sparse(Points, Points); Ax = sparse(Points, Points); else A = zeros(Points); Ax = zeros(Points); end %超过500个点就应用稀疏阵技术

    if getA_Dim == 1
        run1 = 'reshape(no(:,1,:),Points,vad(d));';
    else
        run1 = 'sub2ind(N);'; for k = 1:getA_Dim; run1 = insertBefore(run1, ');', [',reshape(no(:,', num2str(k), ',:),Points,vad(d))']); end
        run2 = 'sub2ind(N);'; for k = 1:getA_Dim; run2 = insertBefore(run2, ');', [',reshape(no(:,', num2str(k), ',:),Points,vad(d)*vad(s))']); end
    end

    for d = 1:getA_Dim

        if drift_yesno(d) %漂移项函数存在，处理它
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            vv = va{d}{1}(ijk(:, d), :);
            no = ijk + permute(reshape(vv(:) * e(d, :), Points, vad(d), getA_Dim), [1, 3, 2]);
            no2 = eval(run1);
            cc = -va{d}{2}(ijk(:, d), :);
            z = (1:Points)' * ones(1, vad(d));
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            Ax = spconvert([z(:), no2(:), cc(:); Points, Points, 0]);
            totf = totf + 1; disp(['getA is running  ', num2str(totf), ' / ', num2str(3 * totfun)]);
            A = A + Ax .* a(:, d).'; % 更新概率密度矩阵
        end

        if diffu_yesno(d, d) %扩散项函数存在，处理它
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

            if diffu_yesno(d, s) %交叉项函数存在，处理它
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

function type = funtype(fun) %根据函数输出数据判断类型
    global getA_domain getA_Dim;
    type = 0;
    x = (0:20)' * diff(getA_domain') / 20 + getA_domain(:, 1)'; %x是[21,getA_Dim]维试错数据
    y = []; y1 = []; y2 = [];

    for no = 1:1
        try [y1, y2] = feval(fun, x); break; end
        try [y1, y2] = feval(fun, x(10, :)); break; end
        try y = feval(fun, x); break; end
        try y = feval(fun, x(10, :)); break; end
    end

    for no = 1:1
        try if ~(isequal(size(y1), [21, getA_Dim]) && isequal(size(y2), [getA_Dim, getA_Dim, 21])); error('*'); end; type = 1; break; end %1型指同时返回[Points,getA_Dim]维和[getA_Dim,getA_Dim,Points]维数据的函数（能同时计算所有漂移系数和所有扩散系数的函数）
        try if ~isequal(size(y), [getA_Dim, getA_Dim, 21]); error('*'); end; type = 2; break; end %2型指返回[getA_Dim,getA_Dim,Points]维数据的函数（能计算所有扩散系数的函数）
        try if ~isequal(size(y), [21, getA_Dim]); error('*'); end; type = 3; break; end %3型指返回[Points,getA_Dim]维数据的函数（能计算所有漂移系数的函数）
        try if ~(isequal(size(y), [21, 1]) | isequal(size(y), [1, 21])); error('*'); end; type = 4; break; end %4型指返回[Points,1]维（或[1,Points]维）数据的函数（能计算漂移系数或计算扩散系数的函数）
        try if ~((isequal(size(y1), [1, getA_Dim]) | isequal(size(y1), [getA_Dim, 1])) && isequal(size(y2), [getA_Dim, getA_Dim])); error('*'); end; type = 5; break; end %5型指同时返回[getA_Dim,1]维（或[1,getA_Dim]维）和[getA_Dim,getA_Dim]维数据的函数（能同时计算所有漂移系数和所有扩散系数的函数）
        try if ~isequal(size(y), [getA_Dim, getA_Dim]); error('*'); end; type = 6; break; end %6型指返回[getA_Dim,getA_Dim]维数据的函数（能计算所有扩散系数的函数）
        try if ~(isequal(size(y), [1, getA_Dim]) | isequal(size(y), [getA_Dim, 1])); error('*'); end; type = 7; break; end %7型指返回[getA_Dim,1]维（或[1,getA_Dim]维）数据的函数（能计算所有漂移系数的函数）
        try if ~isequal(size(y), [1, 1]); error('*'); end; type = 8; break; end %8型指返回[1,1]维数据的函数（能计算1个漂移系数或1个扩散系数的函数）
    end

    if type == 0; error('The input function has no corresponding type.'); end
end
