%���������һ������������(x,y)��������������ǡ�
%�����뵥�����꣬�������Ƕ�����(-pi,pi)֮�䣬�����3����Ϊ(-pi,-pi/2)֮�䣬��1����Ϊ(0,pi/2)֮�䣬��2����Ϊ(pi/2,pi)֮�䣬�ȵ�..
%���ǣ������������꣬���ǰһ�������ڵ�2���ޣ���һ����������3���޵����Σ�����Զ�����2*pi����������������ǡ�
%����x=[1;-1;-1;1]; y=[1;1;-1;-1];
%ph=anglexy(x,y);
%subplot(1,2,1); plot(x,y,'-k'); xlim([-1.2,1.2]); ylim([-1.2,1.2]);
%subplot(1,2,2); plot(1:length(x),ph,'-ok');
%ylim([0,2*pi]); yticks([0 pi/2 pi 3/2*pi 2*pi]); yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
function p=anglexy(x,y)
x=x(:); y=y(:);
p=sign(y).*(1-sign(x))/2*pi+atan(y./x);
pos=find(diff(p)-pi>0);
for k=1:length(pos)
    p(pos(k)+1:end)=p(pos(k)+1:end)-2*pi;
end
pos=find(diff(p)+pi<0);
for k=1:length(pos)
    p(pos(k)+1:end)=p(pos(k)+1:end)+2*pi;
end
end