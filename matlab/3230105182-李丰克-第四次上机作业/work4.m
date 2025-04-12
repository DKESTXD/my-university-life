syms x;
y=exp(-x*x);
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

function X=separate(d,a,b)%均分点
    X=[];
    h=(b-a)/(d-1);
    for i=1:d
        x=a+h*(i-1);
        X=[X,x];
    end    
end

X=separate(100,-5,5);
disp(X);
Y=[];
for i=1:length(X)
    y0=subs(y,x,X(i));
    Y=[Y,y0];
end    
f1=lagrange(x,X,Y);
f2=newton(x,X,Y);
f3=cut(x,X,Y);
x=-5:0.01:5;
f1_cal=matlabFunction(f1);
f2_cal=matlabFunction(f2);
f3_cal=matlabFunction(f3);
f4_cal=matlabFunction(y);
f1_y=f1_cal(x);
f2_y=f2_cal(x);
f3_y=f3_cal(x);
f4_y=f4_cal(x);
figure;
plot(x,f1_y,'r');
hold on;
plot(x,f2_y,'g');
plot(x,f3_y,'b');
plot(x,f4_y,'k');
legend('f1(x)=拉格朗日', 'f2(x)=牛顿', 'f3(x)=分段线性','f4(x)=原函数');
hold off;

figure;
subplot(4,1,1);
plot(x,f1_y,'r');
title('f1(x)=拉格朗日');
subplot(4,1,2);
plot(x,f2_y,'g');
title('f2(x)=牛顿');
subplot(4,1,3);
plot(x,f3_y,'b');
title('f3(x)=分段线性');
subplot(4,1,4);
plot(x,f4_y,'k');
title('f4(x)=原函数');