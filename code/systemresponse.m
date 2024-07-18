%%%%%%%%%%%%%以下程序验证函数systemresponse%%%%%%%%%%
% randn('seed',4321);
% Tstep=0.1; Numtstep=100; Dn1=0.5; Dn2=0.2;
% Ex1=whitenoise(Tstep,Numtstep,Dn1);
% Ex2=whitenoise(Tstep,Numtstep,Dn2);
% initialY=[2,-1];
% sys=@(t,y,w)[0 1;-2 -0.2]*y+[0;w(1)+w(2)*sin(t)];  %系统函数为匿名函数形式
% initial_time=13.58; initial_noise=[0.2145,1.582];
% y1=systemresponse(sys,[Ex1,Ex2],Tstep,initialY,initial_time,initial_noise);
% figure;
% subplot(1,3,1); plot((1:Numtstep)*Tstep+initial_time,y1); title('随机激励下系统响应（systemresponse算得）');
% y2=systemresponse(sys,zeros(Numtstep,2),Tstep,initialY);
% 注意：当系统函数写成非匿名函数形式（如function dydt=sys(t,y,w)....)时，如果放在同一程序下，
% 那么参数sys调用格式为@sys，如果不在同一程序下，而是另为sys.m文件时，那么参数sys调用格式为@sys或"sys"
% subplot(1,3,2); plot((1:Numtstep)*Tstep,y2); title('确定性系统响应（systemresponse算得）');
% %以下通过把噪声赋为零，转成确定性系统，可供ode45使用
% sys1=@(t,y)sys(t,y,[0,0,0]);  %注：其中[0,0,0]即是把噪声赋为零，零值个数大于等于原随机激励系统的噪声数即可;
% subplot(1,3,3); ode45(sys1,[0,Numtstep*Tstep],initialY); title('确定性系统响应（ode45算得）');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time_list, Ys_systemresponse] = systemresponse(varargin)
    % global Ys_systemresponse k_systemresponse
    if (nargin > 6 || nargin < 4); error('Wrong input parameters in systemresponse!'); end
    system_name = varargin{1}; noise = varargin{2}; tstep = varargin{3}; initial_Y = varargin{4};
    initial_time = 0;
    initial_noise=zeros(size(noise,2),1);
    if nargin==6; initial_noise=varargin{6}; end        %如果输入了第6个参数，就作为初始噪声值，否则默认为零
    if nargin >= 5 && isempty(varargin{5}) == 0; initial_time = varargin{5}; end %如果输入了第5个参数，就作为初始时间，否则默认为零
    %system_name是运动系统函数名
    %例：一个受两个相加噪声激励的线性系统，x两点+c*x一点+k*x=noise(1)+noise(2)*sin(t);
    %令y(1)=x，y(2)=x一点，那么运动系统函数格式为：
    %function dydt=sys(t,y,w)     %注：若函数做成了M文件，则调用时要用@sys
    %k=2; c=0.2;
    %dydt=[0 1;-k -c]*y+[0;w(1)+w(2)*sin(t)];
    %end
    %noise是按列排列的噪声序列，即noise(:,1)是第1个噪声，noise(:,2)是第2个噪声，依次类推
    %tstep是计算时间步长
    %initial_Y是初始运动状态，行向量或列向量都行
    %initial_time是初始时刻的时间
    %initial_noise是初始时刻，即initial_time时刻的噪声值，行向量或列向量都行，注意噪声值数量length(initial_noise)与noise列数必须相同
    %输出系统运动响应矩阵Ys的行数与噪声noise行数相同，行向量即是系统状态
    %Ys(i,j)即是initial_time+i*tstep时刻运动状态向量中第j个状态的响应
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    total_tstep = size(noise, 1); %总的计算步长即是噪声序列的长度，即noise的行数，最终运动状态的时间为initial_time+tstep*total_tstep
    Y = initial_Y(:);
    W = initial_noise(:);
    Ys_systemresponse = zeros(total_tstep, length(Y)); %用来贮存系统状态
    total_time = initial_time + tstep*total_tstep;
    time_list = linspace(initial_time, total_time, total_tstep);

    for k_systemresponse = 1:total_tstep
        t = k_systemresponse * tstep;
        W1 = noise(k_systemresponse, :)';
        dk1 = feval(system_name, initial_time + t - tstep, Y, W);
        Y1 = Y + tstep / 2 * dk1;
        dk2 = feval(system_name, initial_time + t - tstep / 2, Y1, (W + W1) / 2);
        Y1 = Y + tstep / 2 * dk2;
        dk3 = feval(system_name, initial_time + t - tstep / 2, Y1, (W + W1) / 2);
        Y1 = Y + tstep * dk3;
        dk4 = feval(system_name, initial_time + t, Y1, W1);
        Y = Y + tstep / 6 * (dk1 + 2 * dk2 + 2 * dk3 + dk4);
        W = W1;
        Ys_systemresponse(k_systemresponse, :) = Y';
        %此处的时刻为initial_time+t，系统状态为Y，列向量;
    end

end

