A=[-2.8  -1.4   0     0;
    1.4   0     0     0;
   -1.8  -0.3  -1.4  -0.6;
     0     0     0.6   0];
B=[1; 0; 1; 0];
C=[0 0 0 1];
D=0;
sys = ss(A, B, C, D);
G=tf(sys)

[z,p,k]=zpkdata(G,'v');
disp('零点:');
disp(z);
disp('极点:');
disp(p);

% 判断稳定性（看极点实部）
if all(real(p) < 0)
    disp('系统稳定');
else
    disp('系统不稳定');
end