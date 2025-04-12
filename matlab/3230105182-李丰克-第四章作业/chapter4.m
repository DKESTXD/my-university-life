syms x;
function f=lagrange(x,X,Y)
    n=length(X);
    f=0;
    for i=1:n
        up=1;
        down=1;
        for j=1:n
            if i~=j
                up=up*(x-X(j));%求基函数分子
            end    
        end
        for j=1:n
            if i~=j
                down=down*(X(i)-X(j));%求基函数分母
            end    
        end
        lk=up/down;%基函数
        f=f+Y(i)*lk;%累加
    end
end

function result=different(X,Y)
    n=length(X);
    result=0;
    for i=1:n
        down=1;
        for j=1:n
            if i~=j
                down=down*(X(i)-X(j));
            end    
        end 
        lk=Y(i)/down;
        result=result+lk;
    end    
end

function f=newton(x,X,Y)
    f=Y(1);
    n=length(X);
    for i=2:n
        k=1;
        result=different(X(1:i),Y(1:i));
        k=result;
        for j=1:i-1
            k=k*(x-X(j));
        end
        f=f+k;
    end
end

function f=cut(x,X,Y)
    f_set=[];
    n=length(X);
    for i=1:n-1
        X_part=X(i:i+1);
        Y_part=Y(i:i+1);
        f_element=lagrange(x,X_part,Y_part);
        f_set=[f_set,f_element];
    end
    f=f_set(1).*(x>=X(1)&x<=X(2));
    for i=2:length(f_set)
        f=f+f_set(i).*(x>X(i)&x<=X(i+1));
    end
end

function f=spline(x,X,Y)
X=X(:);
Y=Y(:);
n=length(X);
h=diff(X);
A=zeros(n,n);
b=zeros(n,1);
A(1,1)=1;
A(n,n)=1;
for i=2:n-1
    A(i,i-1)=h(i-1);
    A(i,i)=2*(h(i-1)+h(i));
    A(i,i+1)=h(i);
    b(i)=3*((Y(i+1)-Y(i))/h(i)-(Y(i)-Y(i-1))/h(i-1));
end
M=A\b;% 求解线性方程组，得到二阶导数M
f=0;
for i=1:n-1 %构建分段符号函数
    segFunc=Y(i)+M(i)/2*(x-X(i))^2+(Y(i+1)-Y(i)-M(i)*h(i)^2/2)/h(i)*(x-X(i))+(M(i+1)-M(i))/(6*h(i))*(x-X(i))^3;
    f=f+segFunc*heaviside(x-X(i))*heaviside(X(i+1)-x);
end
f=f+(Y(n)+M(n)/2*(x-X(n))^2)*heaviside(x-X(n));
end

function result=getadd(m,n,X)
    result=0;
    for i=1:m
        result=result+power(X(i),n);
    end
end

function f=fitting1(x,X,Y)
    a=0;b=0;
    A=zeros(2,2);
    B=zeros(2,1);
    m=length(X);
    A(1,1)=m;
    A(1,2)=getadd(m,1,X);
    A(2,1)=A(1,2);
    A(2,2)=getadd(m,2,X);
    B(1,1)=getadd(m,1,Y);
    result=0;
    for i=1:m
        result=result+X(i)*Y(i);
    end
    B(2,1)=result;
    M=A\B;
    a=M(1);b=M(2);
    f=1/(a+b/x);
end

function f=fitting2(x,X,Y)
    a=0;b=0;
    A=zeros(2,2);
    B=zeros(2,1);
    m=length(X);
    A(1,1)=m;
    A(1,2)=getadd(m,1,X);
    A(2,1)=A(1,2);
    A(2,2)=getadd(m,2,X);
    B(1,1)=getadd(m,1,Y);
    result=0;
    for i=1:m
        result=result+X(i)*Y(i);
    end
    B(2,1)=result;
    M=A\B;
    a=exp(M(1));b=M(2);
    f=a*exp(b/x);
    fprintf("%.16f\n",a);
    fprintf("%.16f\n",b);
end

X=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
Y=[6.7,8.2,9.58,9.5,9.7,10,9.96,9.99,10.49,10.59,10.6,10.8,10.6,10.9,10.76];

X_1=[];
Y_1=[];
X_2=[];
Y_2=[];
for i=1:15
    X_1=[X_1,1/X(i)];
    Y_1=[Y_1,1/Y(i)];
end    
for i=1:15
    X_2=[X_2,1/X(i)];
    Y_2=[Y_2,log(Y(i))];
end    

y1=lagrange(x,X,Y);
y2=newton(x,X,Y);
y3=cut(x,X,Y);
y4=spline(x,X,Y);
y5=fitting1(x,X_1,Y_1);
y6=fitting2(x,X_2,Y_2);
disp(y6);

y1_handle=matlabFunction(y1);
y2_handle=matlabFunction(y2);
y3_handle=matlabFunction(y3);
y4_handle=matlabFunction(y4);
y5_handle=matlabFunction(y5);
y6_handle=matlabFunction(y6);

u=2:0.01:16;
y1_y=y1_handle(u);
y2_y=y2_handle(u);
y3_y=y3_handle(u);
y4_y=y4_handle(u);
y5_y=y5_handle(u);
y6_y=y6_handle(u);

y5_fit=y5_handle(X);
y5_fit_diff=Y-y5_fit;
[a_1,index_1]=max(abs(y5_fit_diff));
y5_fit_max=y5_fit_diff(index_1);
mes_1=mean((Y-y5_fit).^2);
fprintf("第一个拟合方式最大偏差%.16f\n",y5_fit_max);
fprintf("第一个拟合方式均方误差%.16f\n",mes_1);

y6_fit=y6_handle(X);
y6_fit_diff=Y-y6_fit;
[a_2,index_2]=max(abs(y6_fit_diff));
y6_fit_max=y6_fit_diff(index_2);
mes_2=mean((Y-y6_fit).^2);
fprintf("第二个拟合方式最大偏差%.16f\n",y6_fit_max);
fprintf("第二个拟合方式均方误差%.16f\n",mes_2);

figure;
scatter(X,Y);
hold on;
plot(u,y1_y);
plot(u,y2_y);
plot(u,y3_y);
plot(u,y4_y);
legend("散点图", "拉格朗日插值","牛顿插值","分段线性","三次样条插值");
hold off;

figure;
subplot(3,2,1);
scatter(X,Y);
title("散点图");
subplot(3,2,2);
plot(u,y1_y);
title("拉格朗日插值");
subplot(3,2,3);
plot(u,y2_y);
title("牛顿插值");
subplot(3,2,4);
plot(u,y3_y);
title("分段线性");
subplot(3,2,5);
plot(u,y4_y);
title("三次样条插值");

figure;
scatter(X,Y);
hold on;
plot(u,y5_y);
plot(u,y6_y);
legend("散点值","第一种拟合","第二种拟合");
