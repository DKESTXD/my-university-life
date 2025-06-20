num=[1 6 13];
den=conv([1 3 0],[1 1]);
sys=tf(num,den);
rlocus(sys);
K_range=logspace(-2,8,10000);
[P,K]=rlocus(sys,K_range);
P=P';
xi=0.6;
K1=[];
tool=[];
for i=1:length(K)
    for j=1:3
        p=P(i,j);
        xi_tool=-real(p)/abs(p);
        if abs((xi_tool-xi)/xi)<=0.0005
            K1=[K1,K(i)];
            tool=[tool,i];
        end
    end
end
disp(K1);
disp(tool);
k1=(K1(2)+K1(3))/2;
k2=K1(5);
p1=[P(tool(1),1),P(tool(1),2),P(tool(1),3)];
p2=[P(tool(5),1),P(tool(5),2),P(tool(5),3)];
z=[-3+2i,-3-2i,-1];
G1=zpk(z,p1,k1)
G1_real=feedback(tf([k1 6*k1 13*k1],[1 3 0]),tf([1],[1 1]))
G2=zpk(z,p2,k2)
G2_real=feedback(tf([k2 6*k2 13*k2],[1 3 0]),tf([1],[1 1]))