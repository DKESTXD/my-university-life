

![9160943aab45a0e101c865ac53a3f93](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9160943aab45a0e101c865ac53a3f93.png)

![7782313fa188630a0572b1e5ecfc017](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7782313fa188630a0572b1e5ecfc017.png)

# 		 	  						**本科生实验报告**

​								课程名称：   自动控制原理  

​								      姓  名：     李丰克    

​								      学  号：   3230105182   

​                                                                      学  院：  控制科学与工程学院 

​                                                                      专  业：    自动化控制   

​                                                                 指导老师：           

​								 实验日期：    2025.5.7   


























![img](data:image/png;base64,R0lGODlhNwABAHcAMSH+GlNvZnR3YXJlOiBNaWNyb3NvZnQgT2ZmaWNlACH5BAEAAAAALAAAAAA3AAEAgAAAAP///wIIRIynyesN3ysAOw==)![8a4149ef21087a1798b1e3c45e49974](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\8a4149ef21087a1798b1e3c45e49974.png)

 课程名称：<u>自动控制原理</u>       		  指导老师：         		 成    绩：       

 实验名称：<u>MATLAB实验二</u>         	    实验类型：          	       同组学生姓名：       

## 1 实验内容

### 1.1 实验内容1

非单位反馈控制系统的传递函数为：
$$
G(s)=\frac{10A(s^2+8s+20)}{s(s+4)}\ \ \ H(s)=\frac{0.2}{s+2}
$$
绘制系统的根轨迹，确定具有最小阻尼比ξ的放大系数A，并用零、极点、增益形式表示闭环传递函数。

### 1.2 实验内容2

$$
G(s)=\frac{K(s^2+6s+13)}{s(s+3)}\ \ \ H(s)=\frac{1}{s+1}
$$

假设峰值Mp=1.0948，确定满足Mp的ξ值对应的K值，并用零极点增益方式表示闭环传递函数。（计算精度±0.05%）

## 2 代码实现以及结果

### 2.1 实验内容1

开环传递函数为
$$
G(s)H(s)=\frac{2A(s^2+8s+20)}{s(s+4)(s+2)}\\
$$
令2A=k

```matlab
num=[1 8 20];%分子
den=conv([1 4 0],[1 2]);%分母
sys=tf(num,den);
rlocus(sys);%根轨迹
[P,K]=rlocus(sys);%K是增益数；P是每一个增益K对应的极点数，每一个增益对应三个极点
P=P';%由于P原来是3×63的，不方便索引，于是转置
xi=1;
k=0;
tool=0;
for i=1:length(K)
    for j=1:3
        p=P(i,j);
        xi_tool=-real(p)/abs(p);
        if xi_tool<xi%找到最小阻尼比
            xi=xi_tool;
            k=K(i);
            tool=i;
        end
    end
end
disp(xi);
disp(k);
```

根轨迹为

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\decd395c44d19e23dfe006b48378aca.png" alt="decd395c44d19e23dfe006b48378aca" style="zoom:67%;" />

最小阻尼比以及对应增益为

![abe2f9584bc88f90f2d95635071dff6](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\abe2f9584bc88f90f2d95635071dff6.png)

由于K=2A，A=1.0390

利用零点、极点、增益形式表示闭环传递函数。

- 确定零点：由于是非单位反馈，因此闭环传递函数的零点是G(s)的零点与H(s)的极点，为[-4+2i，-4-2i，-2]
- 确定极点：因为特征方程不变，所以此时根轨迹里面对应的极点就是闭环极点，可以用P(tool,i)索引
- 确定增益：因为闭环传递函数的公式为$$\frac{G(s)}{1+G(s)H(s)}$$，那么其增益与G(s)一致，为10A=5k=10.390

```matlab
p_g=[P(tool,1),P(tool,2),P(tool,3)];%闭环极点
z_g=[-4+2i,-4-2i,-2];%闭环零点
G=zpk(z_g,p_g,5*k)
GG=tf([5*k 40*k 100*k],[1 4 0]);%G(s)
HH=tf([0.2],[1 2]);%H(s)
G_1=feedback(GG,HH)
```

此外用定义的方法验证结果，用feedback函数，输出前向通路上的增益函数GG与反馈支路上的增益函数HH，输出为

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\79cbcfcb106b88838cb9397330ffd79.png" alt="79cbcfcb106b88838cb9397330ffd79" style="zoom:67%;" />

两者相同。

### 2.2 实验内容2 

峰值为1.0948，超调量0.0948，$$\sigma=e^{\frac{-\pi \xi}{\sqrt{1-\xi^2}}}$$

代入求解$$\frac{\xi}{\sqrt{1-\xi^2}}=-\frac{ln(0.0948)}{\pi}$$，解得$$\xi=0.6$$

和上一题类似的解法，遍历每一个增益K对应的极点对应的阻尼比，找到满足条件的阻尼比就存下来

```matlab
num=[1 6 13];
den=conv([1 3 0],[1 1]);
sys=tf(num,den);
rlocus(sys);
K_range=logspace(-2,8,10000);%定义高分辨率K的范围
[P,K]=rlocus(sys,K_range);
P=P';
xi=0.6;
K1=[];%存储满足的K
tool=[];%存储对应的下标
for i=1:length(K)
    for j=1:3
        p=P(i,j);%每个增益对应三个极点
        xi_tool=-real(p)/abs(p);
        if abs((xi_tool-xi)/xi)<=0.0005%满足精度百分之五的要求
            K1=[K1,K(i)];
            tool=[tool,i];
        end
    end
end
disp(K1);
disp(tool);
```

由于默认的rlocus函数返回的K只有62个，即步长过大，导致分辨率不够，可能会跳过满足条件的点，因此设置K的范围为(-2,8)取10000个值遍历。

根轨迹为

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\fa2ebb268352f95ad1414c0e6d85b68.png" alt="fa2ebb268352f95ad1414c0e6d85b68" style="zoom:67%;" />

输出为

![af588344de7237f0dfacacd2839a596](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\af588344de7237f0dfacacd2839a596.png)

上面一行为符合的增益K值，下面为对应的下标

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\3f178cbd2a2f71d7d125ca9a0c68662.png" alt="3f178cbd2a2f71d7d125ca9a0c68662" style="zoom:33%;" /><img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\3a3201ea99d56a0dd439caa0ca72ff0.png" alt="3a3201ea99d56a0dd439caa0ca72ff0" style="zoom:33%;" />

从图上找到对应的点，结果大致符合

因为有计算精度0.05%的要求，所以判断条件是一个区间，因此会有多个相近值的出现，实际上只有两个准确值，此处大致认为是K=0.1926与6.2095

```matlab
k1=(K1(2)+K1(3))/2;%K1=0.1926
k2=K1(5);%K2=6.2095
p1=[P(tool(1),1),P(tool(1),2),P(tool(1),3)];%K1对应的极点
p2=[P(tool(5),1),P(tool(5),2),P(tool(5),3)];%K2对应的极点
z=[-3+2i,-3-2i,-1];%零点
G1=zpk(z,p1,k1)
G1_real=feedback(tf([k1 6*k1 13*k1],[1 3 0]),tf([1],[1 1]))
G2=zpk(z,p2,k2)
G2_real=feedback(tf([k2 6*k2 13*k2],[1 3 0]),tf([1],[1 1]))
```

闭环函数与上一题一致，也分别用feedback函数进行了验证，两者大致近似。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c89cf6fbdc6cdc9ce9aa059b24d2b04.png" alt="c89cf6fbdc6cdc9ce9aa059b24d2b04" style="zoom:50%;" /><img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\2adb81b142b1461467725ee5962643b.png" alt="2adb81b142b1461467725ee5962643b" style="zoom:50%;" />

## 3 实验心得体会

通过本实验，学习了matlab中根轨迹的绘制以及其返回值P与K的用法。也学习了用zpk函数求零点、极点、增益方式表示的闭环函数。在求解第一题时，无意将非单位反馈的问题等价成了开环函数为G(s)H(s)的问题，这是**错误的**。从定义角度两者的闭环传递函数分别为$$\frac{G(s)}{1+G(s)H(s)}$$与$$\frac{G(s)H(s)}{1+G(s)H(s)}$$区别在于H(s)，但是在用根轨迹的零点、极点、增益表示闭环传递函数时，两者的区别具体体现在闭环零点与增益上，非单位反馈的闭环零点需要考虑G(s)的零点与H(s)的极点，而单位反馈只需要考虑G(s)的零点，如果忽视了这一点，那么在求闭环传递函数时容易无意中把非单位反馈的问题等价成了开环函数为G(s)H(s)的问题。加深了对非单位反馈的闭环零点的认识。

