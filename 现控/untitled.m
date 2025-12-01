% 绘制传递函数 G(s) = K(s-0.5)/(s-1)^2 的根轨迹
clear;
clc;
close all;

% 定义传递函数 G(s) = K(s-0.5)/(s-1)^2
num = [1, -0.5];   % 分子多项式系数 (s - 0.5)
den = [1, -2, 1];  % 分母多项式系数 (s-1)^2 = s² - 2s + 1

% 绘制根轨迹
figure('Name','根轨迹图','NumberTitle','off');
rlocus(num, den);
grid on;
title('G(s) = K(s-0.5)/(s-1)^2 的根轨迹');
xlabel('实部');
ylabel('虚部');

% 标注极点和零点
hold on;
poles = roots(den);    % 计算极点
zeros = roots(num);    % 计算零点

% 绘制极点（红色圆圈）
plot(real(poles), imag(poles), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
% 绘制零点（绿色叉号）
plot(real(zeros), imag(zeros), 'gx', 'MarkerSize', 10, 'LineWidth', 2);

% 添加图例
legend('根轨迹', '极点', '零点', 'Location', 'Best');

% 计算并标注临界稳定点（穿越虚轴的点）
[K, poles_crit] = rlocfind(num, den);  % 允许用户点击查找特定点的K值
fprintf('临界稳定点的K值: %.4f\n', K);
fprintf('临界稳定点: %.4f + j%.4f\n', real(poles_crit(1)), imag(poles_crit(1)));

% 标注K=0.5时的极点位置（系统刚好稳定）
K_stable = 0.5;
poles_stable = roots(den + K_stable * num);
plot(real(poles_stable), imag(poles_stable), 'bs', 'MarkerSize', 8, 'LineWidth', 2);
text(real(poles_stable(1))+0.1, imag(poles_stable(1))+0.1, 'K=0.5', 'Color', 'b');
