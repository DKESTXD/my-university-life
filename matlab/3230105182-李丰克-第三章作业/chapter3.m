syms x;
function [X,h]=separate(d,a,b)%均分点
    X=[];
    h=(b-a)/(d-1);
    for i=1:d
        x=a+h*(i-1);
        X=[X,x];
    end    
end
[X,h]=separate(11,0,120);
n=length(X);%点数
A=zeros(n-2,n);%定义初始0矩阵
B=[];%常数向量
for i=2:n-1
    k=5*X(i)-X(i)*X(i)/24;
    k=k*h*h/24000000;%常数项
    B=[B,k];
    A(i-1,i-1)=1;
    A(i-1,i)=-2;
    A(i-1,i+1)=1;
end    
A_re=A(:, 2:n-1);
X_re=X(2:n-1);

function x=Gauss(a,b)
    [m,n]=size(a);
    for k=1:n-1
        for i=k+1:n
            factor=a(i,k)/a(k,k);
            for j=k+1:n
                a(i,j)=a(i,j)-factor*a(k,j);
            end
            b(i)=b(i)-factor*b(k);
        end
    end
    x(n)=b(n)/a(n,n);
    for i=n-1:-1:1
        sum=b(i);
       for j=i+1:n
           sum=sum-a(i,j)*x(j);
       end
       x(i)=sum/a(i,i);
    end
end



function x=ColMajorGauss(a,b)
    [m,n]=size(a);
    for k=1:n-1
        col=a(k:n,k);%取出第k列
        col=abs(col);%取绝对值
        [max_value,max_index]=max(col);%获得最大值的索引并转为真实索引
        if max_index+k-1 ~=k%对比是否为当前行，如果不是则换行
            a([k,max_index+k-1], :) = a([max_index+k-1,k],:);
            b([k,max_index+k-1])=b([max_index+k-1,k]);%常数项也要交换
        end
        for i=k+1:n%选好主元后再进行高斯消元
            factor=a(i,k)/a(k,k);
            for j=k:n
                a(i,j)=a(i,j)-factor*a(k,j);
            end
            b(i)=b(i)-factor*b(k);
        end
    end
    x(n)=b(n)/a(n,n);
    for i=n-1:-1:1
        sum=b(i);
       for j=i+1:n
           sum=sum-a(i,j)*x(j);
       end
       x(i)=sum/a(i,i);
    end
end

function x=cholesky(A,b)
    n=length(b);L=zeros(n,n);y=zeros(n,1);x=zeros(n,1);
    for i=1:n%分解出L
        for k=1:i
            sumLk2=0;
            if k==i
                for j=1:k-1
                    sumLk2=sumLk2+L(k,j)^2;
                end
                L(k,k)=sqrt(A(k,k)-sumLk2);
            else
                for j=1:k-1
                    sumLk2=sumLk2+L(i,j)*L(k,j);
                end
                L(i,k)=(A(i,k)-sumLk2)/L(k,k);
            end
        end
    end
    for i=1:n %前代求解
        sumLy=0;
        for j=1:i-1
            sumLy=sumLy+L(i,j)*y(j);
        end
        y(i)=(b(i)-sumLy)/L(i,i);
    end
    for i=n:-1:1%后代求解
        sumLx=0;
        for j=i+1:n
            sumLx=sumLx+L(j,i)*x(j);
        end
        x(i)=(y(i)-sumLx)/L(i,i);
    end
    x=x';
end

function x=thomas(A, f)
    n=size(A,1);
    c=zeros(1,n-1); % 上对角线元素
    b=zeros(1,n);     % 主对角线元素
    a=zeros(1,n-1); % 下对角线元素
    b(:)=diag(A);
    for i = 1:n-1
        c(i) = A(i+1, i);
    end
    for i = 2:n
        a(i-1) = A(i, i-1);
    end
    be=zeros(1,n-1);
    be(1)=c(1)/b(1);
    for i=2:n-1
        be(i)=c(i)/(b(i)-a(i-1)*be(i-1));
    end
    y=zeros(1,n);
    y(1)=f(1)/b(1);
    for i=2:n
        y(i)=(f(i)-a(i-1)*y(i-1))/(b(i)-a(i-1)*be(i-1));
    end
    x=zeros(1,n);
    x(n)=y(n);
    for i=n-1:-1:1
        x(i)=y(i)-be(i)*x(i+1);
    end    
end

function [x,k]=jacobi(a,b,x0,epsilon)
    n=size(a,1);
    x=zeros(1,n);
    flag=0;
    k=0;
    while 1
    k=k+1;
    for i=1:n
        sum=0;
        for j=1:n
            if i~=j
                sum=sum+a(i,j)*x0(j);
            end
        end
        x(i)=(b(i)-sum)/a(i,i);
        if abs((x(i)-x0(i))/x(i))<epsilon
            flag=1;
        end    
        x0(i)=x(i);
    end    
    if flag==1
        break;
    end    
    end
end

function [x,k]=GS(a,b,x0,epsilon)
    n=size(a,1);
    x=zeros(1,n);
    flag=0;
    k=0;
    while 1
    k=k+1;
    for i=1:n
        sum1=0;
        sum2=0;
        if i~=1
            for j=1:i-1
                sum1=sum1+a(i,j)*x(j);
            end
        end
        for j=i+1:n
            sum2=sum2+a(i,j)*x0(j);
        end
        x(i)=(b(i)-sum1-sum2)/a(i,i);
        if abs((x(i)-x0(i))/x(i))<epsilon
            flag=1;
        end    
        x0(i)=x(i);
    end    
    if flag==1
        break;
    end   
end
end

x0=[1,1,1,1,1,1,1,1,1];
y1=Gauss(A_re,B);
y2=ColMajorGauss(A_re,B);
y3=cholesky(A_re,B);
y4=thomas(-A_re,-B);
[y5,k5]=jacobi(A_re,B,x0,eps);
[y6,k6]=GS(A_re,B,x0,eps);


u=0:0.01:120;
y_theroy=120*x*x*x/(12*12*24000000)-x*x*x*x/(12*24*24000000)-x*120*120*120/(12*24*24000000);
y1=[0,y1,0];
y2=[0,y2,0];
y3=[0,y3,0];
y4=[0,y4,0];
y5=[0,y5,0];
y6=[0,y6,0];
y_theroy_handle=matlabFunction(y_theroy);
y_theroy_y=y_theroy_handle(u);

figure;
subplot(3, 2, 1);
plot(X,y1);
title("原始高斯消元");
subplot(3, 2, 2);
plot(X,y2);
title("列主元高斯消元");
subplot(3, 2, 3);
plot(X,y3);
title("Cholesty分解");
subplot(3, 2, 4);
plot(X,y4);
title("Thomas算法");
subplot(3, 2, 5);
plot(X,y5);
title("雅各比迭代");
subplot(3, 2, 6);
plot(X,y6);
title("GS迭代");

figure;
plot(X, y1, 'r', 'DisplayName', "原始高斯消元");
hold on; % 保持当前图形，允许在同一图上绘制多个函数
plot(X, y2, 'g', 'DisplayName', "列主元高斯消元");
plot(X, y3, 'b', 'DisplayName', "Cholesty分解");
plot(X, y4, 'c', 'DisplayName', "Thomas算法");
plot(X, y5, 'y', 'DisplayName', "雅各比迭代");
plot(X, y6, 'm', 'DisplayName', "GS迭代");
plot(u, y_theroy_y, 'k', 'DisplayName', '理论解');
legend show;
hold off;

function show(y)
    for i=1:length(y)
        fprintf("%.16f\n",y(i));
    end 
    fprintf("\n");
end
fprintf("原始高斯\n");
show(y1);
fprintf("列主元\n");
show(y2);
fprintf("cholesty分解\n");
show(y3);
fprintf("Thomas算法\n");
show(y4);
fprintf("雅各比迭代\n");
show(y5);
fprintf("GS迭代\n");
show(y6);
fprintf("雅各比迭代\n");
disp(k5);
fprintf("GS迭代\n");
disp(k6);