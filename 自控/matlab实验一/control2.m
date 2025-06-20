t=0:0.01:10;
figure;
hold on;
for xi=0.1:0.1:2.0
    wn=6;
    sys=wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    [y,t]=step(sys, t);
    plot(t,y,'DisplayName',['\xi=' num2str(xi)]);
end
title('单位阶跃响应(omega_n=6)');
legend show;
grid on;
hold off;

figure;
hold on;
wn_values=[2, 4, 6, 8, 10, 12];
for wn=wn_values
    xi=0.7;
    sys=wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    [y,t]=step(sys, t);
    plot(t,y,'DisplayName',['\omega_n =' num2str(wn)]);
end
title('单位阶跃响应(xi=0.7)');
legend show;
grid on;
hold off;
