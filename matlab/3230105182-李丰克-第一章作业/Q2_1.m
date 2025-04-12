%思路一
function result=answer1(x)
    result=(1-cos(x))/sin(x);
end
%思路二
function result=answer2(x)
    sum_sin=0;
    sum_cos=0;
    term=0;
    i=0;
    while 1
        term=((-1)^i)*(x^(2*i+1))/factorial(2*i+1);
        sum_sin=sum_sin+term;
        i=i+1;
        if abs(term)<eps
            break;
        end
    end
    fprintf("%.16f\n",sum_sin);
    i=0;
    term=0;
    while 1
        term=((-1)^i)*(x^(2*i))/factorial(2*i);
        sum_cos=sum_cos+term;
        i=i+1;
        if abs(term)<eps
            break;
        end
    end
    fprintf("%.16f\n",sum_cos);
    result=(1-sum_cos)/sum_sin;
end
%思路三
function result=answer3(x)
    result=tan(x/2);
end
%思路三点一
function result=answer4(x)
    result=x+(x^3)/3+(x^5)*2/15+(x^7)*7/315+(x^9)*62/2835+(x^11)*1382/155925+(x^13)*21844/6081075+(x^15)*929569/638512875
end
x1=-0.01;
x1_single=single(x1);
x1_double=double(x1);
x2=0.0001;
x2_single=single(x2);
x2_double=double(x2);

y1_ans1_single=answer1(x1_single);
y1_ans1_double=answer1(x1_double);
y2_ans1_single=answer1(x2_single);
y2_ans1_double=answer1(x2_double);
fprintf("%.16f\n",y1_ans1_single);
fprintf("%.16f\n",y1_ans1_double);
fprintf("%.16f\n",y2_ans1_single);
fprintf("%.16f\n",y2_ans1_double);
fprintf("\n");

y1_ans2_single=answer2(x1_single);
y1_ans2_double=answer2(x1_double);
y2_ans2_single=answer2(x2_single);
y2_ans2_double=answer2(x2_double);
fprintf("%.16f\n",y1_ans2_single);
fprintf("%.16f\n",y1_ans2_double);
fprintf("%.16f\n",y2_ans2_single);
fprintf("%.16f\n",y2_ans2_double);
fprintf("\n");

y1_ans3_single=answer3(x1_single);
y1_ans3_double=answer3(x1_double);
y2_ans3_single=answer3(x2_single);
y2_ans3_double=answer3(x2_double);
fprintf("%.24f\n",y1_ans3_single);
fprintf("%.24f\n",y1_ans3_double);
fprintf("%.24f\n",y2_ans3_single);
fprintf("%.24f\n",y2_ans3_double);
fprintf("\n");

y1_ans4_single=answer4(x1_single/2);
y1_ans4_double=answer4(x1_double/2);
y2_ans4_single=answer4(x2_single/2);
y2_ans4_double=answer4(x2_double/2);
fprintf("%.24f\n",y1_ans4_single);
fprintf("%.24f\n",y1_ans4_double);
fprintf("%.24f\n",y2_ans4_single);
fprintf("%.24f\n",y2_ans4_double);
fprintf("\n");