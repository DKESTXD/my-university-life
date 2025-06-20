

![d7cc8bdafcbb1c353f7196b1813f6c1](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\d7cc8bdafcbb1c353f7196b1813f6c1.png)

![660d70829dd553575a8613d06abac43](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\660d70829dd553575a8613d06abac43.png)

# 		 	  		    		        **本科生实验报告**

​								课程名称：   自动控制原理  

​								      姓  名：     李丰克    

​								      学  号：   3230105182   

​                                                                      学  院：  控制科学与工程学院 

​                                                                      专  业：    自动化控制   

​                                                                 指导老师：           

​								 实验日期：    2025.6.4   























![1d19afe2be866f433710e5af5229af2](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\1d19afe2be866f433710e5af5229af2.png)

 课程名称：<u>自动控制原理</u>       		  指导老师：         		 成    绩：       

 实验名称：<u>MATLAB实验三</u>         	    实验类型：          	       同组学生姓名：       

## 1 实验内容

### 1.1 实验内容1

单位反馈开环系统$$H(s)=\frac{50}{(s+1)(s+5)(s-2)}$$，绘制系统Nyquist曲线，判断系统闭环稳定性，绘制出闭环系统的脉冲响应

### 1.2 实验内容2

控制系统的传递函数为$$G(s)=\frac{K(s+1)}{(s+2)(s^2+4s+5)}$$，用对数频率特性确定相位裕度大于45°时的Km值

## 2 代码实现以及结果

### 2.1 实验内容1

```matlab
num=[50];
den=[1 4 -7 -10];
sys=tf(num,den)
nyquist(sys);
```

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20250605132729346.png" alt="image-20250605132729346" style="zoom:67%;" />

由于是顺时针N=-1，且在右侧极点数为1，PR=1，$$Z_R=P_R-N=2$$，因此不稳定

```matlab
G=feedback(sys,1);
figure;
impulse(G);
```

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\6cf8afc319525b2b0912fcd580e19fa.png" alt="6cf8afc319525b2b0912fcd580e19fa" style="zoom:60%;" />

### 2.2 实验内容2 

 

```matlab
den=conv([1 2],[1 4 5]);
K_values=[];
K_scan=linspace(1, 100, 100000);
for i=1:length(K_scan)
    K_val=K_scan(i);
    sys_temp=tf([K_val K_val],den);
    % 获取Bode图的相位裕度
    [Gm,Pm,Wcg,Wcp]=margin(sys_temp);
    if Pm>45
        K_values=[K_values, K_val];
    end
end

if ~isempty(K_values)
    fprintf('%.4f\n',max(K_values));
end
```

输出值为43.6877，故满足相位裕度大于45°的最大Km值为43.6877

```
num=[43.6877 43.6877];
den=conv([1 2],[1 4 5]);
sys=tf(num,den);
[Gm,Pm,Wcg,Wcp]=margin(sys);
disp(Pm);
```

验证一下

当K=43.6877时：![e67babceba6fd8530a45cca40858fba](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\e67babceba6fd8530a45cca40858fba.png)

当K=43.6878时：![c593223de332c4efaa06861c286719e](D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c593223de332c4efaa06861c286719e.png)

## 3 实验心得体会

通过本次实验，学习了如何用matlab画系统的nyquist图象，加深了对nyquist稳定性判据的印象。此外也学习了用matlab求解相位裕度、幅值裕度的方法，已经如何利用其求解系统参数。
