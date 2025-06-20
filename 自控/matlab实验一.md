

![9160943aab45a0e101c865ac53a3f93](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9160943aab45a0e101c865ac53a3f93.png)

![7782313fa188630a0572b1e5ecfc017](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7782313fa188630a0572b1e5ecfc017.png)

# 		 	  **本科生实验报告**

​								课程名称：   自动控制原理  

​								      姓  名：     李丰克    

​								      学  号：   3230105182   

​                                                                      学  院：  控制科学与工程学院 

​                                                                      专  业：    自动化控制   

​                                                                 指导老师：           

​								 实验日期：    2025.4.9   













 
![img](data:image/png;base64,R0lGODlhNwABAHcAMSH+GlNvZnR3YXJlOiBNaWNyb3NvZnQgT2ZmaWNlACH5BAEAAAAALAAAAAA3AAEAgAAAAP///wIIRIynyesN3ysAOw==)![8a4149ef21087a1798b1e3c45e49974](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\8a4149ef21087a1798b1e3c45e49974.png)

 课程名称：<u>自动控制原理</u>       		  指导老师：         		 成    绩：       

 实验名称：<u>MATLAB实验一</u>         	    实验类型：          	       同组学生姓名：       

## 1 实验内容

### 1.1 实验内容1

<img src="D:\新建文件夹 (5)\img1\微信截图_20250409165229.png" alt="微信截图_20250409165229" style="zoom:50%;" />

方块图如图，求输入输出传递函数。（并与方框图法得到的传递函数进行比较）

### 1.2 实验内容2

1，典型二阶系统$$H(s)=\frac{\omega_n^2}{s^2+2\zeta\omega_n+\omega_n^2}$$，其中$$\omega_n$$是自然频率(无阻尼振荡频率)，ξ为相对阻尼系数，试绘制：
	1）当$$\omega_n$$=6，ξ分别为0.1,0.2,0.3,...,2.0时的单位阶跃响应(绘制在一张图上)

​	2）当ξ=0.7，$$\omega_n$$取2,4,6,8,10,12时的单位阶跃响应(绘制在一张图上)

2，编程计算二阶系统$$G(s)=\frac{1}{s^2+s+1}$$的时域指标（上升时间，超调量，峰值时间，稳态时间）

## 2 代码实现以及结果

### 2.1 实验内容1

要用matlab模拟方块图，要有Control System Toolbox

画出信号流图
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\ce8fc57be007bc8fb4fb7323543c0d3.jpg" alt="ce8fc57be007bc8fb4fb7323543c0d3" style="zoom: 33%;" />

对每条支路的传递函数做定义(包括负号)并且编号。

```matlab
s=tf("s");
G1=1/(s+1);
G2=1/(s+2);
G3=1/(s+3);
G4=-1/(s+4);
G5=1/(s+5);
G6=-1;
sys=append(G1,G2,G3,G4,G5,G6);
Q=[1 5 6;
    2 1 4;
    3 2 0;
    4 3 0;
    5 2 0;
    6 3 0];
inputs=1;
outputs=3;
sysc=connect(sys,Q,inputs,outputs);
G=tf(sysc)
```

Q矩阵反映了各支路的连接情况，每一行对应每一条支路，第一位为该支路的输出编号，后面若干位为该支路的输入。

第一行：支路1的输入为5，6的输出；
第二行：支路2的输入为1，4的输出；
第三行：支路3的输入为2的输出；
第四行：支路4的输入为3的输出；
第五行：支路5的输入为2的输出；
第六行：支路6的输入为3的输出；

输出为：
![微信截图_20250410211727](D:\新建文件夹 (5)\img1\微信截图_20250410211727.png)

即$$G(s)=\frac{s^2+9s+20}{s^5+15s^4+85s^3+226s^2+282s+133}$$

用方框图法解得的传递函数与matlab解得的一致
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7cb791ac2c9aa540f379bba0538be03.jpg" alt="7cb791ac2c9aa540f379bba0538be03" style="zoom:67%;" />

### 2.2 实验内容2 

```matlab
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
```

输出为
<img src="D:\新建文件夹 (5)\img1\微信截图_20250409170143.png" alt="微信截图_20250409170143" style="zoom:50%;" /><img src="D:\新建文件夹 (5)\img1\微信截图_20250409170155.png" alt="微信截图_20250409170155" style="zoom:50%;" />

```matlab
s=tf('s');
G=1/(s^2+s+1);
[y,t]=step(G);
info=stepinfo(G);
fprintf("上升时间%.8f\n",info.RiseTime);
fprintf("超调量%.8f%%\n",info.Overshoot);
fprintf("峰值时间%.8f\n",info.PeakTime);
fprintf("稳态时间%.8f\n",info.SettlingTime);
figure;
plot(t, y);
title('单位阶跃响应');
grid on;

```

输出为：

<img src="D:\新建文件夹 (5)\img1\微信截图_20250409170641.png" alt="微信截图_20250409170641" style="zoom:67%;" />![aa9400866861850a627c67e7cd4d09a](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\aa9400866861850a627c67e7cd4d09a.png)

## 3 实验心得体会

control system toolbox库涵盖了课程中大部分的内容，如方框图、信号流图、时域响应指标等等，其中计算传递函数相当方便，只是定义方框图比较麻烦，尤其是Q矩阵，要先画出信号流图，再编号，再每条支路写出输入输出，得到矩阵。