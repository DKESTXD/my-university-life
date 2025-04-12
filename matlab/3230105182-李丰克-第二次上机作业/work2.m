function [k,p,err,P]=fixpt(g,p0,tol,max1)
%input-g是迭代函数，p0是不动点初始值，tol是容限，max1是最大迭代次数
%output-k是当前迭代的次数，err是误差，p是不动点最终结果，P是p迭代的序列
P(1)=p0;
for k=2:max1
	P(k)=feval(g,P(k-1));
	err=abs(P(k)-P(k-1));
	relerr=err/abs(P(k)+eps);
	p=P(k);
	if (err<tol)||(relerr<tol),break;end
end
end

function show(k,p,err,P)
disp(['迭代次数: ', num2str(k)]);
fprintf("值%.16f\n",p);
fprintf("误差上界%.28f\n",err);
fprintf("\n");
for i=1:k
    fprintf("%.16f\n",P(i));
end
fprintf("\n");
end

f=@(x)(power(x,4)-5*x-10);
g1=@(x)((power(x,4)-10)/5);
g2=@(x)(power((10+5*x),0.25));
g3=@(x)(nthroot((10+5*x)./x,3));
g4=@(x)(sqrt((10+5*x)./(x.^2)));

function P_value=cycle(P,g)
Pnew=P(1); 
for i=2:length(P)
    Pnew=[Pnew, P(i), P(i)]; 
end
P_value(1)=feval(g,Pnew(1));
for i=2:length(Pnew)
    if mod(i,2)==0
        P_value(i)=P_value(i-1);
    end
    if mod(i,2)==1
        P_value(i)=feval(g,P_value(i-1));
    end
end    
    plot(Pnew,P_value,"-o");
    hold on;
    x=[-2:0.01:-0.01,0.01:0.01:2.5];
    g_value=g(x);
    tool=@(x)(x);
    tool_value=tool(x);
    plot(x,g_value);
    plot(x,tool_value);
    hold off
end


figure;
x=-5:0.01:5;
y_value=f(x);
plot(x,y_value);

p0=2;
p0=double(p0);

[k1,p1,err1,P1]=fixpt(g1,p0,eps,1000);
[k2,p2,err2,P2]=fixpt(g2,p0,eps,1000);
[k3,p3,err3,P3]=fixpt(g3,p0,eps,1000);
[k4,p4,err4,P4]=fixpt(g4,p0,eps,1000);


figure;
subplot(4,1,1);
x1=1:length(P1);
plot(x1,P1,"-o");

subplot(4,1,2);
x2=1:length(P2);
plot(x2,P2,"-o");

subplot(4,1,3);
x3=1:length(P3);
plot(x3,P3,"-o");

subplot(4,1,4);
x4=1:length(P4);
plot(x4,P4,"-o");

figure
cycle(P1,g1); 

figure
cycle(P2,g2); 

figure
cycle(P3,g3); 

figure
cycle(P4,g4); 

show(k1,p1,err1,P1);
show(k2,p2,err2,P2);
show(k3,p3,err3,P3);
show(k4,p4,err4,P4);
disp(f(p1));
disp(f(p2));