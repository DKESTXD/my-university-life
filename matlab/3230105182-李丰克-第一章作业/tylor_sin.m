x=0.5;
x_single=single(x);
x_double=double(x);
eps_four=1e-5;
eps_single=eps("single");
eps_double=eps("double");

function sum=sin_tylor(x,eps)
    i_value=[];
    term_value=[];
    sum_value=[];
    sum=0;
    term=0;
    i=0;
    while 1
        term=((-1)^i)*(x^(2*i+1))/factorial(2*i+1);
        sum=sum+term;
        i_value=[i_value;i];
        term_value=[term_value;term];
        sum_value=[sum_value;sum];
        i=i+1;
        if abs(term)<eps
            break;
        end
    end
    T = table(i_value,term_value,sum_value,'VariableNames',{'Index','Term','Sum'});
    disp(T);
end

y_double_four=sin_tylor(x_double,eps_four);
y_single_four=sin_tylor(x_single,eps_four);
y_double=sin_tylor(x_double,eps_double);
y_single=sin_tylor(x_single,eps_single);
fprintf("%.16f\n",y_double_four);
fprintf("%.16f\n",y_single_four);
fprintf("%.16f\n",y_double);
fprintf("%.16f\n",y_single);

y_double_four=vpa(y_double_four,4);
y_single_four=vpa(y_single_four,4);
y_double=vpa(y_double,4);
y_single=vpa(y_single,4);

disp(y_double_four);
disp(y_single_four);
disp(y_double);
disp(y_single);

y_single_double=sin_tylor(x_double,eps_single);
y_double_single=sin_tylor(x_double,eps_single);
fprintf("%.16f\n",y_single_double);
fprintf("%.16f\n",y_double_single);