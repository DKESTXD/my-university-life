double x;
double x_new;
x=input("输入值");
a=x;
while 1
    x_new=(2*x+a/(x^2))/3;
    if abs(x_new-x)<eps %eps
        break;
    end
    disp(abs(x_new-x));
    x=x_new;
end
disp(abs(x_new-x));
x=vpa(x_new,5);
disp(x);