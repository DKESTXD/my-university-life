num=[1 8 20];
den=conv([1 4 0],[1 2]);
sys=tf(num,den);
rlocus(sys);
[P,K]=rlocus(sys);
P=P';
xi=1;
k=0;
tool=0;
for i=1:length(K)
    for j=1:3
        p=P(i,j);
        xi_tool=-real(p)/abs(p);
        if xi_tool<xi
            xi=xi_tool;
            k=K(i);
            tool=i;
        end
    end
end
disp(xi);
disp(k);
p_g=[P(tool,1),P(tool,2),P(tool,3)];
z_g=[-4+2i,-4-2i,-2];
G=zpk(z_g,p_g,5*k)
GG=tf([5*k 40*k 100*k],[1 4 0]);
HH=tf([0.2],[1 2]);
G_1=feedback(GG,HH)
