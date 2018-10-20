
x=-2*pi:2*pi/40:2*pi;
q=ones(length(x),1)*4*pi/length(x);
q=q/sum(q);

subplot(3,2,1)
y1=0.06+sin(x)*0.01;
p1=y1/sum(y1);
dkl1=sum((p1+eps).*log((p1./q')+eps));
plot(x,y1,'k','Linewidth',2)
legend(strcat('DKL=',num2str(dkl1)))
ylim([0 0.1])
%xlim([-3*pi/2 3*pi/2])

subplot(3,2,2)
y2=0.06+sin(x)*0.015;
p2=y2/sum(y2);
dkl2=sum((p2+eps).*log((p2./q')+eps));
plot(x,y2,'k','linewidth',2)
ylim([0 0.1])
legend(strcat('DKL=',num2str(dkl2)))
%xlim([-3*pi/2 3*pi/2])

subplot(3,2,3)
y3=0.06+sin(x)*0.01+0.01*rand(size(x));
p3=y3/sum(y3);
dkl3=sum((p3+eps).*log((p3./q')+eps));
plot(x,y3,'k','Linewidth',2)
ylim([0 0.1])
legend(strcat('DKL=',num2str(dkl3)))
%xlim([-3*pi/2 3*pi/2])

subplot(3,2,4)
y4=0.06+sin(x)*0.01+0.02*rand(size(x));
p4=y4/sum(y4);
dkl4=sum((p4+eps).*log((p4./q')+eps));
plot(x,y4,'k','Linewidth',2)
ylim([0 0.1])
legend(strcat('DKL=',num2str(dkl4)))
%xlim([-3*pi/2 3*pi/2])

subplot(3,2,5)
y5=0.01*sawtooth(x)+0.06;
p5=y5/sum(y5);
dkl5=sum((p5+eps).*log((p5./q')+eps));
plot(x,y5,'k','Linewidth',2)
ylim([0 0.1])
legend(strcat('DKL=',num2str(dkl5)))
%xlim([-3*pi/2 3*pi/2])


subplot(3,2,6)
y6=0.01*sawtooth(x)+0.06+0.01*rand(size(x));
p6=y6/sum(y6);
dkl6=sum((p6+eps).*log((p6./q')+eps));
plot(x,y6,'k','Linewidth',2)
ylim([0 0.1])
legend(strcat('DKL=',num2str(dkl5)))
%xlim([-3*pi/2 3*pi/2])