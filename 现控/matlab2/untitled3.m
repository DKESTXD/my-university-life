A=[-2 -3; 4 9];
B=[3; 1];
% 指定期望闭环极点
p = [-1+2j, -1-2j];

% 方法一：acker
K1 = acker(A, B, p)

% 方法二：place（更稳定）
K2 = place(A, B, p)
Acl=A-B*K1
% 验证闭环极点
eig(A - B*K1)