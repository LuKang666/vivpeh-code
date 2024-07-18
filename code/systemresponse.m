%%%%%%%%%%%%%���³�����֤����systemresponse%%%%%%%%%%
% randn('seed',4321);
% Tstep=0.1; Numtstep=100; Dn1=0.5; Dn2=0.2;
% Ex1=whitenoise(Tstep,Numtstep,Dn1);
% Ex2=whitenoise(Tstep,Numtstep,Dn2);
% initialY=[2,-1];
% sys=@(t,y,w)[0 1;-2 -0.2]*y+[0;w(1)+w(2)*sin(t)];  %ϵͳ����Ϊ����������ʽ
% initial_time=13.58; initial_noise=[0.2145,1.582];
% y1=systemresponse(sys,[Ex1,Ex2],Tstep,initialY,initial_time,initial_noise);
% figure;
% subplot(1,3,1); plot((1:Numtstep)*Tstep+initial_time,y1); title('���������ϵͳ��Ӧ��systemresponse��ã�');
% y2=systemresponse(sys,zeros(Numtstep,2),Tstep,initialY);
% ע�⣺��ϵͳ����д�ɷ�����������ʽ����function dydt=sys(t,y,w)....)ʱ���������ͬһ�����£�
% ��ô����sys���ø�ʽΪ@sys���������ͬһ�����£�������Ϊsys.m�ļ�ʱ����ô����sys���ø�ʽΪ@sys��"sys"
% subplot(1,3,2); plot((1:Numtstep)*Tstep,y2); title('ȷ����ϵͳ��Ӧ��systemresponse��ã�');
% %����ͨ����������Ϊ�㣬ת��ȷ����ϵͳ���ɹ�ode45ʹ��
% sys1=@(t,y)sys(t,y,[0,0,0]);  %ע������[0,0,0]���ǰ�������Ϊ�㣬��ֵ�������ڵ���ԭ�������ϵͳ������������;
% subplot(1,3,3); ode45(sys1,[0,Numtstep*Tstep],initialY); title('ȷ����ϵͳ��Ӧ��ode45��ã�');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time_list, Ys_systemresponse] = systemresponse(varargin)
    % global Ys_systemresponse k_systemresponse
    if (nargin > 6 || nargin < 4); error('Wrong input parameters in systemresponse!'); end
    system_name = varargin{1}; noise = varargin{2}; tstep = varargin{3}; initial_Y = varargin{4};
    initial_time = 0;
    initial_noise=zeros(size(noise,2),1);
    if nargin==6; initial_noise=varargin{6}; end        %��������˵�6������������Ϊ��ʼ����ֵ������Ĭ��Ϊ��
    if nargin >= 5 && isempty(varargin{5}) == 0; initial_time = varargin{5}; end %��������˵�5������������Ϊ��ʼʱ�䣬����Ĭ��Ϊ��
    %system_name���˶�ϵͳ������
    %����һ�������������������������ϵͳ��x����+c*xһ��+k*x=noise(1)+noise(2)*sin(t);
    %��y(1)=x��y(2)=xһ�㣬��ô�˶�ϵͳ������ʽΪ��
    %function dydt=sys(t,y,w)     %ע��������������M�ļ��������ʱҪ��@sys
    %k=2; c=0.2;
    %dydt=[0 1;-k -c]*y+[0;w(1)+w(2)*sin(t)];
    %end
    %noise�ǰ������е��������У���noise(:,1)�ǵ�1��������noise(:,2)�ǵ�2����������������
    %tstep�Ǽ���ʱ�䲽��
    %initial_Y�ǳ�ʼ�˶�״̬��������������������
    %initial_time�ǳ�ʼʱ�̵�ʱ��
    %initial_noise�ǳ�ʼʱ�̣���initial_timeʱ�̵�����ֵ�������������������У�ע������ֵ����length(initial_noise)��noise����������ͬ
    %���ϵͳ�˶���Ӧ����Ys������������noise������ͬ������������ϵͳ״̬
    %Ys(i,j)����initial_time+i*tstepʱ���˶�״̬�����е�j��״̬����Ӧ
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    total_tstep = size(noise, 1); %�ܵļ��㲽�������������еĳ��ȣ���noise�������������˶�״̬��ʱ��Ϊinitial_time+tstep*total_tstep
    Y = initial_Y(:);
    W = initial_noise(:);
    Ys_systemresponse = zeros(total_tstep, length(Y)); %��������ϵͳ״̬
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
        %�˴���ʱ��Ϊinitial_time+t��ϵͳ״̬ΪY��������;
    end

end

