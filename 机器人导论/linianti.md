## 23-24

一，名词解释

机器人：人造的，模仿人类或者其他自然生命体的，能够实现感知，运动，操作，编程等功能的，包含体能和智能两部分组成的机器

传感器：用于定量感知环境里特定物质属性的机械、电子、化学设备，并且能够将各种物理量、化学量精准转化为电信号，再经过电子电路或者计算机进行分析与处理，从而实现对这些量的检测。

二，

2， 控制器结构？冯诺依曼？

输入是从检测线的传感器得到的信息，进入运算器进行运算，得到轮子电机的命令，再通过控制器下发给输出机构。

3， 1/20K*0.4

4， 速度传感器：增量式光电编码器，单位时间数脉冲或者相邻脉冲的时间
	加速度传感器：F=ma
	光纤陀螺仪：方位角传感器
	激光雷达测距传感器：时飞法

5，均匀取点，删去在障碍上的点，在邻域内连线，删去经过障碍的点

6，RRT RRT*

7，VRB  VO   flocking models

VO：感知他机的速度与位置，对每个机器人都计算VOA|B（vb）(解释一下这个公式)，用线性规划得到一个安全速度，执行速度更新位置

8（1）惰性齿轮，等比传动，中间过渡，改变传动方向

（2）d=mz

（3）求传动比的 120

2*（-1）\*（-1）

（4）转速是120rmp  400HZ

9，内参
|fx    s     cx|  
|0     fy     cy|
|0      0      1 |

外参

【R  |   t】

内参包括焦距fx，fy，光学中心cx，cy，像素倾斜系数s

外参包括旋转矩阵和平移矩阵，表示相机在世界坐标系的位姿

内参不变外参变

## 22-23

1， 体能：能量的产生运输与消耗，依赖电荷驱动和传动维持机器的运转

智能：机器人感知，依赖传感器和电子线路实现控制反馈过程

2，冯诺依曼

3，用于定量感知周围环境特定物质属性的电子机械化学设备，并且能够将各种化学量物理量精确转化为电信号，再经过电子电路或计算机分析与处理，实现对这些量的检测

4，循迹小车

```c
const int leftsensor=2;
const int rightsensor=3;
const int leftmotor=9;
const int rightmotor=10;
const highspeed=255;
const lowspeed=150;
void setup{
    pinMode=(leftsensor,INPUT);
    pinMode=(rightsensor,INPUT);
    pinMode=(leftmotor,OUTPUT);
    pinMode=(rightmotor,OUTPUT);
    
}
void loop{
    int leftsensorvalue=digitalRead(leftsensor);
    int rightsensorvalue=digitalRead(rightsensor);
    
    if(leftsensorvalue==HIGH&&rightsensorvalue==LOW){
        analogWrite(leftmotor,lowspeed);//用pwm波控制电机转速
        analogWrite(rightmotor,highspeed);
    }//后面以此类推
}
```

5，Ea=Ke  n

T=Km  I

P=T omega

1/20K*0.6

6，（1）二级减速

（2）20mm

（3）首/末=-3

7

8，采样，PRM，RRT

搜索：dijkstra  A*