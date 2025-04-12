syms x;
f=1/(1+x^2);
N1=100;
N2=10000;
N1_single=single(N1);
N1_double=double(N1);
N2_single=single(N2);
N2_double=double(N2);
function result=answer(N)
    result=atan(N+1)-atan(N);
end

N1_ans1_single=int(f,x,N1_single,N1_single+1);
N1_ans1_double=int(f,x,N1_double,N1_double+1);
N2_ans1_single=int(f,x,N2_single,N2_single+1);
N2_ans1_double=int(f,x,N2_double,N2_double+1);
fprintf("%.24f\n",N1_ans1_single);
fprintf("%.24f\n",N1_ans1_double);
fprintf("%.24f\n",N2_ans1_single);
fprintf("%.24f\n",N2_ans1_double);
fprintf("\n")

N1_ans2_single=answer(N1_single);
N1_ans2_double=answer(N1_double);
N2_ans2_single=answer(N2_single);
N2_ans2_double=answer(N2_double);
fprintf("%.24f\n",N1_ans2_single);
fprintf("%.24f\n",N1_ans2_double);
fprintf("%.24f\n",N2_ans2_single);
fprintf("%.24f\n",N2_ans2_double);