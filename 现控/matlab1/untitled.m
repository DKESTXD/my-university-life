%定义离散系统
num=[0.632];              
den=[1, -1.368, 0.568];
Ts=1;%采样周期1
sys=tf(num,den,Ts);
N=500;%总采样点数
k=0:N-1;
t=k*Ts; 
%生成离散方波输入
T=100;%方波周期10
u=square(2*pi*k/T); %方波输入
% 计算离散系统响应
y=lsim(sys,u,t);

figure;
hold on; grid on;
% 绘制输入信号（离散点）
stem(t,u,'b' );
% 绘制输出响应（离散点）
stem(t,y,'r');
% 设置坐标轴和图例
xlabel('kT (采样周期)');
ylabel('幅值');
legend('输入 u(kT)', '输出 y(kT)');
hold off;
