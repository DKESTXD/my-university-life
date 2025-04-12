u=0:0.01:1;
f_y=humps(u);
plot(u,f_y);
syms x;
f=1/((x-0.3)^2+0.01)+1/((x-0.9)^2+0.04)-6;
f_gauss=(1/((0.5*x+0.2)^2+0.01)+1/((0.5*x-0.4)^2+0.04)-6)/2;
f_handle=matlabFunction(f);
f_gauss_handle=matlabFunction(f_gauss);
df_true=-140;
intf_true=29.858325395498675089500892382438;
function result=error(y_true,y_pred)
    result=abs(y_true-y_pred)/abs(y_true);
end

function result=diff(f,x0,h)
    result=(-f(x0+2*h)+4*f(x0+h)-3*f(x0))/(2*h);
end

function result=diff_three(f,x0,h)
	result=(f(x0+h)-f(x0-h))/(2*h);
end

function result=trape(f,X)
    result=(X(2)-X(1))*(f(X(2))+f(X(1)))/2;
end
function result=simpson(f,X)
    h=(X(3)-X(1))/2;
    result=h*(f(X(1))+4*f(X(2))+f(X(3)))/3;
end
function result=simpson8(f,X)
    h=(X(4)-X(1))/3;
    result=3*h*(f(X(1))+3*f(X(2))+3*f(X(3))+f(X(4)))/8;
end

function X=separate(d,a,b)%均分点
    X=[];
    h=(b-a)/(d-1);
    for i=1:d
        x=a+h*(i-1);
        X=[X,x];
    end    
end

X=separate(121,0,1);
sum1=0;
sum2=0;
sum3=0;
sum4=0;
for(i=1:40)
   sum1=sum1+simpson8(f_handle,X(3*i-2:3*i+1)); 
end
for(i=1:60)
   sum2=sum2+simpson(f_handle,X(2*i-1:2*i+1));
end
for(i=1:120)
   sum3=sum3+trape(f_handle,X(i:i+1));
end
for(i=1:20)
   sum4=sum4+simpson8(f_handle,X(3*i-2:3*i+1));
end
for(i=30:59)
   sum4=sum4+simpson(f_handle,X(2*i+1:2*i+3));
end    

function [T,S,C,R]=romberg(k,f,a,b)
    T(1)=(b-a)*(f(b)+f(a))/2;
    for i=1:k
        sum=0;
        for j=1:2^(i-1)
            sum=sum+f(a+(2*j-1)*(b-a)/(2^i));
        end
        T(i+1)=T(i)/2+(b-a)*sum/(2^i);
    end    
    for i=1:k
        S(i)=(4*T(i+1)-T(i))/3;
    end    
    for i=1:k-1
        C(i)=(16*S(i+1)-S(i))/15;
    end  
    for i=1:k-2
        R(i)=(64*C(i+1)-C(i))/63;
    end  
end

function result=gauss2(f)
    result=f(0.5773502692)+f(-0.5773502692);
end
function result=gauss3(f)
    result=0.5555555556*f(0.7745966692)+0.5555555556*f(-0.7745966692)+0.8888888889*f(0);
end
function result=gauss4(f)
    result=0.3478548451*(f(0.8611363116)+f(-0.8611363116))+0.6521451549*(f(0.3399810436)+f(-0.3399810436));
end
function result=gauss5(f)
    result=0.2369268851*(f(0.9061798459)+f(-0.9061798459))+0.4786286705*(f(0.5384693101)+f(-0.5384693101))+0.5688888889*f(0);
end
function result=gauss6(f)
    result=0.1713244924*(f(0.9324695142)+f(-0.9324695142))+0.3607615730*(f(0.6612093865)+f(-0.6612093865))+0.4679139346*(f(0.2386191861)+f(-0.2386191861));
end


a1=diff(f_handle,0.5,0.001);
a2=diff_three(f_handle,0.5,0.001);
disp("高精度微分公式：")
disp(a1);
disp(error(df_true,a1));
disp("一阶微分三点公式：")
disp(a2);
disp(error(df_true,a2));


disp("二十个辛普森3/8与三十个辛普森")
disp(sum4);
disp(error(intf_true,sum4));
k=4;
[T,S,C,R]=romberg(k,f_handle,0,1);
disp(T);
disp(S);
disp(C);
disp(R);
disp(error(intf_true,T(k+1)));
disp(error(intf_true,S(k)));
disp(error(intf_true,C(k-1)));
disp(error(intf_true,R(k-2)));

result=gauss6(f_gauss_handle);
disp("n=6时")
disp(result);
disp(error(intf_true,result));