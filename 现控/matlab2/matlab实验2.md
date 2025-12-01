

![ScreenShot_2025-10-16_201650_973](D:\新建文件夹 (5)\1\ScreenShot_2025-10-16_201650_973.png)

![ScreenShot_2025-10-16_201700_565](D:\新建文件夹 (5)\1\ScreenShot_2025-10-16_201700_565.png)

# 		 	  						                                                                                                                                                                                                                      **本科生实验报告**

​								课程名称：   现代控制原理  

​								      姓  名：     李丰克    

​								      学  号：   3230105182   

​                                                                      学  院：  控制科学与工程学院 

​                                                                      专  业：    自动化控制   

​                                                                 指导老师：           

​								 实验日期：    2025.10.30  














![ScreenShot_2025-10-16_201753_717](D:\新建文件夹 (5)\1\ScreenShot_2025-10-16_201753_717.png)

 课程名称：<u>现代控制原理</u>       		  指导老师：         		  成    绩：       

 实验名称：<u>MATLAB实验二</u>                     实验类型：                          同组学生姓名：

## 1 实验内容

### 1.1 实验内容1

给定连续系统状态空间方程，求传递函数模型和零极点模型，并判断其稳定性。

![微信图片_20251028125628_386_49](D:\新建文件夹 (5)\1\微信图片_20251028125628_386_49.png)

### 1.2 实验内容2

#### 1.2.1 Question1

已知受控系统，设计状态反馈阵 K，使系统闭环极点为[-1+j2，-1-j2]， (分别采用上课所讲方法直接编程和 matlab 函数 place 或 acker 方法) 

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_130422_376.png" alt="ScreenShot_2025-10-28_130422_376"  />

#### 1.2.2 Question2

已知系统，设计全维状态观测器，使其极点为-3，-4，-5(分别采用书上方法直接编程和 matlab 函数 estim 方法) 

![ScreenShot_2025-10-28_130700_344](D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_130700_344.png)

## 2 代码实现以及结果

### 2.1 实验内容一

```matlab
%建立状态空间模型并转为传递函数
A=[-2.8  -1.4   0     0;
    1.4   0     0     0;
   -1.8  -0.3  -1.4  -0.6;
     0     0     0.6   0];
B=[1; 0; 1; 0];
C=[0 0 0 1];
D=0;
sys = ss(A, B, C, D);
G=tf(sys)

%得到零极点模型
[z,p,k]=zpkdata(G,'v');
disp('零点:');
disp(z);
disp('极点:');
disp(p);

%判断稳定性（看极点实部）
if all(real(p) < 0)
    disp('系统稳定');
else
    disp('系统不稳定');
end
```

输出为：

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_132513_318.png" alt="ScreenShot_2025-10-28_132513_318" style="zoom: 50%;" /><img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_132550_056.png" alt="ScreenShot_2025-10-28_132550_056" style="zoom: 50%;" />

### 2.2 实验内容二

#### 2.2.1 Question1

##### 2.2.1.1 用直接法求

先判断能控性：
$$
rankQ_C=rank\begin{pmatrix} b & Ab  \end{pmatrix} =rank\begin{pmatrix} 3 & -9 \\ 1 & 21 \end{pmatrix}=2
$$
系统能控。

假设闭环Acl=A+BK，用理想闭环特征多项式和实际上的闭环特征多项式det(λI-A-BK)对应系数相等可解k1，k2

```matlab
syms k1 k2 s;
A=[-2 -3; 4 9];
B=[3; 1];
K=[k1,k2];
real=simplify(det(s*eye(2)-A-B*K));%实际闭环多项式
disp('闭环特征多项式为：')
pretty(expand(real));
desire=s^2+2*s+5;%目标闭环多项式
%系数相等求解
eqs=coeffs(expand(real-desire), s);
eqs=fliplr(eqs); %让最高次在前
sol=solve(eqs == 0, [k1 k2]);
K_sol=[sol.k1,sol.k2]

```

运行结果得到

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_134607_811.png" alt="ScreenShot_2025-10-28_134607_811" style="zoom:75%;" />          <img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_134613_878.png" alt="ScreenShot_2025-10-28_134613_878" style="zoom:75%;" />
$$
x^·=(A+BK)x+Bu=\begin{pmatrix} -6.7917 &-15.6250  \\ 2.4028 & 4.7917 \end{pmatrix}x+\begin{pmatrix} 3  \\ 1  \end{pmatrix}u
$$


##### 2.2.1.2 用matlab内置函数

```matlab
A=[-2 -3; 4 9];
B=[3; 1];
%指定期望闭环极点
p=[-1+2j, -1-2j];
%用acker方法
K1=acker(A,B,p)
%用place方法
K2=place(A,B,p)
%验证闭环极点
eig(A - B*K1)
```

反馈阵：						闭环极点：

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_135033_692.png" alt="ScreenShot_2025-10-28_135033_692" style="zoom:75%;" />                     <img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-28_135038_380.png" alt="ScreenShot_2025-10-28_135038_380" style="zoom:75%;" />
$$
x^·=(A-BK)x+Bu=\begin{pmatrix} -6.7917 &-15.6250  \\ 2.4028 & 4.7917 \end{pmatrix}x+\begin{pmatrix} 3  \\ 1  \end{pmatrix}u
$$


两种方法求出的K与直接法求出的数值是一样的，只是符号相反，这是由于这两种方法默认闭环状态方程Acl=A-BK，而刚才的方法是定义的Acl=A+BK，于是会有符号上的不同，但是最终闭环函数的形式是一样的。

#### 2.2.2 Question2

##### 2.2.2.1 用直接法

首先判断能观性

```matlab
A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
Qo=[C; C*A; C*(A^2)];
disp(Qo);
r=rank(Qo);
disp(r);
```

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-30_142721_367.png" alt="ScreenShot_2025-10-30_142721_367" style="zoom:67%;" />

结果能观

```MATLAB
syms l1 l2 l3 s;
A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
L=[l1;l2;l3];
real=simplify(det(s*eye(3)-A+L*C));%实际闭环多项式
pretty(expand(real));
desire=simplify((s+3)*(s+4)*(s+5));%目标闭环多项式
pretty(expand(desire));
%系数相等求解
eqs=coeffs(expand(real-desire), s);
eqs=fliplr(eqs); %让最高次在前
sol=solve(eqs == 0, [l1 l2 l3]);
L_sol=[sol.l1;sol.l2;sol.l3]
```

结果得到

<img src="D:\新建文件夹 (5)\1\微信图片_20251030143625_390_49.png" alt="微信图片_20251030143625_390_49" style="zoom:67%;" />
$$
x^{·}=(A-LC)x+Bu+Ly=\begin{pmatrix} -119 &-120 &0 \\ 103 & 105 & 1\\-210 &-210&2  \end{pmatrix}x+\begin{pmatrix} 120  \\ -103 \\210 \end{pmatrix}y+\begin{pmatrix} 1  \\ 0 \\1 \end{pmatrix}u
$$

##### 2.2.2.2 用matlab内置函数

```matlab
A=[1 0 0; 0 2 1; 0 0 2];
C=[1 1 0];
p=[-3 -4 -5];

%用acker计算L
L=acker(A', C', p)';
%用estim构造估计器
sys=ss(A, [], C, []);
G=estim(sys, L, 1)

disp('L ='); disp(L);
disp('eig(A-L*C) ='); disp(eig(A - L*C));

```

estim用于根据已算出来的L构建完整的状态观测器系统，可用能控行的acker方法来求L，由于能观性与能控性的对偶关系，应该用转置的A与C。

运行结果的L为：			以及状态观测器为：

<img src="D:\新建文件夹 (5)\1\ScreenShot_2025-10-30_145553_414.png" alt="ScreenShot_2025-10-30_145553_414" style="zoom:67%;" />                                         <img src="D:\新建文件夹 (5)\1\微信图片_20251030145648_391_49.png" alt="微信图片_20251030145648_391_49" style="zoom:50%;" />
$$
x^{·}=(A-LC)x+Bu+Ly=\begin{pmatrix} -119 &-120 &0 \\ 103 & 105 & 1\\-210 &-210&2  \end{pmatrix}x+\begin{pmatrix} 120  \\ -103 \\210 \end{pmatrix}y+\begin{pmatrix} 1  \\ 0 \\1 \end{pmatrix}u
$$

## 3 实验心得

​	通过本实验习得了如何用matlab实现传递函数模型和状态空间模型以及零极点模型之间的互相转化并判断稳定性。另外也加深了对状态反馈设计器与状态观测器设计方法的理解，其中状态观测器的求解与状态反馈可共用一个函数，因为对偶原理，输入改为AT和cT即可。
