%本程序根据一段连续的坐标(x,y)，生成连续的相角。
%若输入单个坐标，则算出相角定义在(-pi,pi)之间，比如第3象限为(-pi,-pi/2)之间，第1象限为(0,pi/2)之间，第2象限为(pi/2,pi)之间，等等..
%但是，输入连续坐标，如果前一个坐标在第2象限，下一个坐标进入第3象限的情形，相角自动增加2*pi，以生成连续的相角。
%例：x=[1;-1;-1;1]; y=[1;1;-1;-1];
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