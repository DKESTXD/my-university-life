A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
p=[-3 -4 -5];

% 用 acker 计算 L（注意转置的技巧）
L=acker(A', C', p)';

% 用 estim 构造估计器（或直接验证 A-LC 的特征值）
sys=ss(A, [], C, []);
G=estim(sys, L, 1)    % 这里 estim 需要 Control System Toolbox

disp('L ='); disp(L);
disp('eig(A-L*C) ='); disp(eig(A - L*C));
