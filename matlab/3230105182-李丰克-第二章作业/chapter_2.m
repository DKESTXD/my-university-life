x=0:0.01:10;
f=@(x)(1000000-100000*x-250000*sqrt(x));
y1_value=f(x);
plot(x, y1_value);
grid();

function [c,err,yc,k]=bisect(f,a,b,delta)
%input-f是函数，a、b是左右边界值，delta是容限
%output-c是零点，yc是f(c)，err是误差
ya=feval(f,a);
yb=feval(f,b);
if(ya*yb>0)
	return;
end
max1=1+round((log(b-a)-log(delta))/log(2)); 
for k=1:max1
	c=(a+b)/2;
	yc=feval(f,c);
	if yc==0
		a=c;
		b=c;
	elseif ya*yc<0
		b=c;
		yb=yc;
	else
		a=c;
		ya=yc;
	end
	if b-a<delta
		break;
    end
end
c=(a+b)/2;
err=abs(b-a)/2;
yc=feval(f,c);
end

function [c,err,yc,k]=regular(f,a,b,delta,epsilon,max1)
%input-f是函数，a、b是左右边界点，delta是相邻两次迭代值之间误差的容限，epsilon是当前点的函数值与0之间差的容限，max1是最 %  大迭代次数
%output-c是零点，yc是f(c)，err是误差
ya=feval(f,a);
yb=feval(f,b);
if ya*yb>0
	return;
end
for k=1:max1
	dx=yb*(b-a)/(yb-ya);
	c=b-dx;
	ac=c-a;
	yc=feval(f,c);
	if yc==0 
		break;
	elseif yb*yc>0
		b=c;
		yb=yc;
	else
		a=c;
		ya=yc;
	end
	dx=min(abs(dx),ac);
	if abs(dx)<delta
		break;
	end
	if abs(yc)<epsilon
		break;
	end
end 
dx=yb*(b-a)/(yb-ya);
c=b-dx;
err=abs(b-a);
yc=feval(f,c);
end

function [c,err,yc,k]=regularnew(f,a,b,delta,epsilon,max1)
%input-f是函数，a、b是左右边界点，delta是相邻两次迭代值之间误差的容限，epsilon是当前点的函数值与0之间差的容限，max1是最 %  大迭代次数
%output-c是零点，yc是f(c)，err是误差
ya=feval(f,a);
yb=feval(f,b);
count_a=0;
count_b=0;
if ya*yb>0
	return;
end
for k=1:max1
	dx=yb*(b-a)/(yb-ya);
	c=b-dx;
	ac=c-a;
	yc=feval(f,c);
	if yc==0 
		break;
	elseif yb*yc>0
        count_a=count_a+1;
        count_b=0;
		b=c;
		yb=yc;
    else
        count_b=count_b+1;
        count_a=0;
		a=c;
		ya=yc;
	end
	dx=min(abs(dx),ac);
    if count_a>=2
        ya=ya/2;
        count_a=0;
    end
    if count_b>=2
        yb=yb/2;
        count_b=0;
    end
	if abs(dx)<delta
		break;
	end
	if abs(yc)<epsilon
		break;
	end
end 
err=abs(b-a);
yc=feval(f,c);
end

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

function [p0,err,k,y]=newton(f,df,p0,delta,epsilon,max1)
%input-f是函数，df是函数一阶导，p0是初始值，delta迭代两次差的容限，epsilon是迭代当前结果函数值与0之间差的容限
%max1最大迭代次数
%output-p0是最后结果，k是当前迭代次数，y是f(p0)函数值
for k=1:max1
	p1=p0-feval(f,p0)/feval(df,p0);
	err=abs(p1-p0);
	relerr=2*err/(abs(p1)+delta);
	p0=p1;
	y=feval(f,p0);
	if(err<delta)||(relerr<delta)||(abs(y)<epsilon),break;end
end
end

function [p1,err,k,y]=secant(f,p0,p1,delta,epsilon,max1)
%input-f是函数，p0是初始值，p1也是估计值，delta是相邻两项的容限，epsilon是y的容限，max1是最大迭代次数
%output-p1是估计值，err是误差，y是函数值
for k=1:max1
	p2=p1-feval(f,p1)*(p1-p0)/(feval(f,p1)-feval(f,p0));
	err=abs(p2-p1);
	relerr=2*err/(abs(p2)+delta);
	p0=p1;
	p1=p2;
	y=feval(f,p1);
	if (err<delta)||(relerr<delta)||(abs(y)<epsilon),break;end
end
end

function [p1,err,k,y]=secantnew(f,p0,varepsilon,delta,epsilon,max1)
%input-f是函数，p0是初始值，varepsilon是扰动小量的倍数，delta是相邻两项的容限，epsilon是y的容限
%max1是最大迭代次数
%output-p1是估计值，err是误差，y是函数值，k是迭代次数
for k=1:max1
	p1=p0-varepsilon*p0*feval(f,p0)/(feval(f,p0+varepsilon*p0)-feval(f,p0));
	err=abs(p1-p0);
	relerr=2*err/(abs(p1)+delta);
	p0=p1;
	y=feval(f,p0);
	if(err<delta)||(relerr<delta)||(abs(y)<epsilon),break;end
end
end

xl=4;
xu=5;
xl=double(xl);
xu=double(xu);

g=@(x)(10-2.5*sqrt(x));
x0=5;
x0=double(x0);

df=@(x)(-100000-125000/sqrt(x));

[c1,err1,yc1,k1]=bisect(f,xl,xu,10*eps);
fprintf("二分法结果%.16f\n",c1);
fprintf("二分法误差上界%.28f\n",err1);
fprintf("%.28f\n",yc1);
disp(k1);

[c2,err2,yc2,k2]=regular(f,xl,xu,10*eps,10*eps,1000);
fprintf("试位法结果%.16f\n",c2);
fprintf("试位法误差上界%.28f\n",err2);
fprintf("%.28f\n",yc2);
disp(k2);

[c3,err3,yc3,k3]=regularnew(f,xl,xu,10*eps,10*eps,1000);
fprintf("试位法修正版结果%.16f\n",c3);
fprintf("试位法修正版误差上界%.28f\n",err3);
fprintf("%.28f\n",yc3);
disp(k3);

[k4, c4, err4, P] = fixpt(g, x0, 10*eps, 1000);
fprintf("不动点迭代法结果%.16f\n",c4);
fprintf("不动点迭代法误差上界%.28f\n",err4);
fprintf("%.28f\n",f(c4));
disp(k4);

[c5,err5,k5,yc5]=newton(f,df,x0,10*eps,10*eps,1000);
fprintf("N-R法结果%.16f\n",c5);
fprintf("N-R法结果误差上界%.28f\n",err5);
fprintf("%.28f\n",yc5);
disp(k5);

[c6,err6,k6,yc6]=secant(f,xl,xu,10*eps,10*eps,1000);
fprintf("割线法结果%.16f\n",c6);
fprintf("割线法结果误差上界%.28f\n",err6);
fprintf("%.28f\n",yc6);
disp(k6);

[c7,err7,k7,yc7]=secantnew(f,x0,0.01,10*eps,10*eps,1000);
fprintf("改进割线法结果%.16f\n",c7);
fprintf("改进割线法结果误差上界%.28f\n",err7);
fprintf("%.28f\n",yc7);
disp(k7);

fprintf("%.16f\n",fzero(f,4));