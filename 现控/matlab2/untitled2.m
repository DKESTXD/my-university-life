syms k1 k2 s;
A=[-2 -3; 4 9];
B=[3; 1];
K=[k1,k2];
real=simplify(det(s*eye(2)-A-B*K));%实际闭环多项式
disp('闭环特征多项式为：')
pretty(expand(real));
desire=s^2+2*s+5;%目标闭环多项式
%系数相等求解
eqs=coeffs(expand(real-desire), s);
eqs=fliplr(eqs); %让最高次在前
sol=solve(eqs == 0, [k1 k2]);
K_sol=[sol.k1,sol.k2]
