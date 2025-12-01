A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
Qo=[C; C*A; C*(A^2)];
disp(Qo);
r=rank(Qo);
disp(r);
syms l1 l2 l3 s;
A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
L=[l1;l2;l3];
real=simplify(det(s*eye(3)-A+L*C));%实际闭环多项式
pretty(expand(real));
desire=simplify((s+3)*(s+4)*(s+5));%目标闭环多项式
pretty(expand(desire));
%系数相等求解
eqs=coeffs(expand(real-desire), s);
eqs=fliplr(eqs); %让最高次在前
sol=solve(eqs == 0, [l1 l2 l3]);
L_sol=[sol.l1;sol.l2;sol.l3]