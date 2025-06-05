## 1 绪论

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\1e53cf2b4512d44853985eb18ed2705.jpg" alt="1e53cf2b4512d44853985eb18ed2705" style="zoom: 15%;" />

### 1.1 分类

**工业机器人**：焊接机器人、磨抛加工机器人、激光加工机器人、喷涂机器人、搬运机器人、机床机器人、冲压机器人、真空机器人等
**服务机器人**：扫地机器人、拖地机器人、擦窗机器人、陪伴型机器人、教育机器人、休闲娱乐机器人、送餐机器人、迎宾机器人、酒店机器人、商场导购机器人、银行柜台机器人、巡检机器人、下水道工作机器人、深海工作机器人、微型机器人、教育机器人、室内安保机器人、室外巡逻机器人、汽车/飞机清洗机器人、消防机器人、管路勘探机器人、导游机器人、公共场所清洁服务机器人、地雷探测机器人、无人驾驶机器人、太空探测及其热、反恐防暴机器人、小型侦查机器人。

**非工业机器人就是服务机器人**，自动坦克也是服务机器人，因为国防用途的机器人也属于服务机器人

### 1.2 定义

机器人的**定义**：人造的，模仿人类或者自然生命体的机器，能够有实现感知、操作、运动、编程等的能力。

本质：人造的机器，有人类的特性（体能，智能）

- 体能：能量存储，消耗，传输。依赖电荷、驱动和传动 维持机器的运转
- 智能：机器人感知，中央处理器，信息传输与处理。依赖传感器和电子线路实现控制反馈过程

核心技术

- 机器人的肌肉——驱动
- 机器人的骨骼——机构
- 机器人的运动——建模及控制
- 机器人的感官——传感器
- 机器人的知觉——识别理解
- 机器人的作业——决策规划
- 机器人的协作——多智能体

## 2 嵌入式

### 2.1 冯诺依曼

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\3deb392f1257601b010534ea2e62f18.png" alt="3deb392f1257601b010534ea2e62f18" style="zoom:67%;" />

审题：
有哪几个模块：输入设备，输出设备，运算器，存储器，控制器
or要求画出结构图：注意分数据流和控制流即可。

按照冯诺依曼结构具体分析某个机器人

### 2.2 编程题

Arduino编程：结构：

```c
#define led 13
void setup{
	pinMode(led,OUTPUT);
}
void loop{
	digitalWrite(led,HIGH);
	delay(1000);//延迟1s,单位为毫秒
	digitalWrite(led,LOW);
	delay(1000);
}
```

此为LED闪烁的程序，arduino有setup和loop两个主要函数。

```c
void setup(){
	Serial.begin(9600);//设置波特率9600
}
void loop(){
	Serial.println("jiqirendaolun manji");//输出内容
	delay(1000);
}
```

如何实现LED灯不闪烁但是亮度只有正常的20%。

```c
#define led 13
Void setup（）
{
	pinMode(led,OUTPUT); //设定led管脚为输出引脚
}
Void loop()
{
	digitalWrite（led,HIGH); //设置led为高电平，点亮led
	delay(2);//延迟2ms
    digitalWrite（led,LOW); //设置led为低电平，熄灭led
	delay(8);//延迟8ms
}
```

由于人眼刷新率有限，所以只要延迟时间够短就可以。可以推广到输出电压调节(电容滤波)，直流电机控制，信号传输(舵机控制)。
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\918a259230099760f1629fa2b005d8d.png" alt="918a259230099760f1629fa2b005d8d" style="zoom:50%;" />

PWM波

```c
void setup()
{
    pinMode(13,OUTPUT);//设定13号端口为输出
}

void loop()
{
    digitalWrite(13,HIGH);
    delayMicroseconds(100);//大约10%占空比的1KHZ方波
    digitalWrite(13,LOW);
    delayMicroseconds(900);
}
```

- delay()/delayMicroseconds()：用于延时，第一个单位为毫秒，第二个为微秒。  

- analogWrite():模拟 I/O 口输出，一般用于 PWM 输出，如：

  analogWrite(13,127)，为在13 号引脚处输出一个占空比为 50%的 PWM 方波，后一参数 0 表示关， 255 表示全开  

0——占空比0%
64——占空比25%
127——占空比50%
191——占空比75%
255——占空比100%

## 3 传感器

### 3.1 定义

用于**定量**感知环境**特定物质属性**的**电子、机械、化学**设备，并能够把各种物理量和化学量等**精确**地变换为**电信号**，再经由电子电路或计算机进行**分析与处理**，从而对这些量进行检测.

### 3.2 分类(举例)

**内部传感器**：测量机器人自身状态，常用于底层运动控制。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\49845a4968f1c8dd6b1d106c74ff683.png" alt="49845a4968f1c8dd6b1d106c74ff683" style="zoom:50%;" />

**外部传感器**：测量机器人所处环境，部分用于底层运动控制，部分用于上层运动规划

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7258826c426c0a2260b6ab7927f5ccc.png" alt="7258826c426c0a2260b6ab7927f5ccc" style="zoom:50%;" />

### 3.3 特性

**静态特性**：指检测系统的输入为**不随时间变化**的恒定信号时，系统的输出与输入之间的关系

基本要求：输出相对于输入保持一定的对应关系

- 灵敏度 ：对输入信号变化的响应敏感度（**越高越好**）
- 信噪比（ S/N ）：传感器输出信号中信号分量与噪声分量的平方平均值之比
- 线性 ：输入与输出量之间为线性比例关系
- 时滞（回差）：输入量滞后，输出量也滞后
- 稳定性： 输入量恒定，输出量向一个方向偏移（温漂/零漂）
- 精度：
  - **准确度**：测量值对真值的偏离程度；
  - **精密度**：测量相同对象，每次测量会得到不同测量值



**动态特性**：指检测系统的输入为随时间变化的信号时，系统的输出与输入之间的关系  

- 瞬态响应特性
- 频率响应特性  

**传感器选择**：

- 测量条件
- 传感器的性能
- 成本、尺寸、重量
- 使用条件

### 3.4 主要传感器及其基本原理

- 运动传感器 
- 方位角传感器 
- 力觉传感器 
- 接触觉传感器 
- 接近觉传感器 
- 定位和测距传感器

#### 3.4.1 运动传感器 

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\aded4154a850bf6575c2197aa65084d.png" alt="aded4154a850bf6575c2197aa65084d" style="zoom:50%;" />

##### 3.4.1.1 电位器

原理就是滑动变阻器。

旋转式：测量角位移 
直线式：测量线位移

- 单独使用

- 和其他传感器（如编码器）一起使用

  - 用电位器检测起始位置

  - 用编码器检测关节和连杆的当前位置

##### 3.4.1.2 编码器

根据测量介质分：光电码盘，磁编码器
根据测量结果分：增量式，绝对式

**增量式光电编码器**：

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7d8ed5b15cd80a8e1e6e976c547b0f2.png" alt="7d8ed5b15cd80a8e1e6e976c547b0f2" style="zoom: 33%;" />

有两个码道，等间距透明的栅格，有光源持续照射栅格，在旋转过程中，根据**数脉冲数可以得到转过的角度**。
顺逆时针方向，根据两条码道谁先提前来判断。
	设内部码道为A，外部码道为B，**逆时针外部比内部领先90°，顺时针外部比内部落后90°**。

**绝对式光电编码器**：

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\0ece77ed96bf1476ae29010633f7273.png" alt="0ece77ed96bf1476ae29010633f7273" style="zoom:50%;" />

每一个位置对应一个二进制编码，只要转到那里，就输出那个位置的绝对编码。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\eb17c9620556de75c4ab8236a0f21a0.png" alt="eb17c9620556de75c4ab8236a0f21a0" style="zoom:50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\6eb644c2263da9711cf747c9d20e6a4.png" alt="6eb644c2263da9711cf747c9d20e6a4" style="zoom:50%;" />

**磁编码器**：
非接触，适合用在连续旋转物体上。

##### 3.4.1.3 速度传感器

用编码器测量速度：统计指定时间内脉冲信号数量。采样时间越短，越解决瞬时速度。

当速度快时，计频率，即指定时间内经过的脉冲数量。
当速度慢时，计周期，即相邻脉冲的时间间隔。

##### 3.4.1.4 加速度传感器

基本原理：利用加速度造成某个介质产生变形，通过测量其变形量并用相关电路转化成电压输出，如压电晶体

F=ma

#### 3.4.2 方位角传感器

用于测量机器人的方向和倾角，可进行机器人位姿估计。

主要传感器：指南针，陀螺仪，倾角仪。

##### 3.4.2.1 指南针

- 原理：地球南北极的磁力
- 缺点：易受其他磁性物质和人类环境的干扰

##### 3.4.2.2 陀螺仪
**机械陀螺仪**：

- 原理：运用物体高速旋转时，角动量很大，旋转轴会一直稳定指向一个方向的性质
- 缺点：工艺要求很高，结构复杂，精度受到多方面的制约

**光纤陀螺仪**：

- 原理：光束的速度保持不变，光干涉原理
- 优点：结构紧凑，灵敏度高，工作可靠
- 缺点：价格高

**MEMS陀螺仪**：

- 原理：地球自转运动而作用于地球上运动质点的偏向力
- 优点：体积小、重量轻；低成本；低功耗；大量程；
- 缺点：**yaw角不准**

**惯导**（important）陀螺仪

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c79f5179d6d049601be6b1ef7556bf9.png" alt="c79f5179d6d049601be6b1ef7556bf9" style="zoom: 33%;" />

#### 3.4.3 力觉传感器

压阻：半导体压阻效应
压电：压电效应：某些晶体电压与压力有关系
压容：电容机理：电容量由电极面积和两个电极间的距离决定，压力改变间距

##### 力矩传感器

原理：当力矩作用在弹性轴上，轴会产生扭曲变形，存在剪切应变和应力

#### 3.4.4 测距传感器

##### 3.4.4.1 超声波接近觉传感器
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\07c0f51544ec70c8273a1c5c110cd82.png" alt="07c0f51544ec70c8273a1c5c110cd82" style="zoom:50%;" />

**原理**：利用压电传感器生成声波，采用测量传输时间法测距

- **优点**：
  - 穿透性强
  - 波比较慢，时间长，处理方便


- **缺点**：
  - 声波传输速度低，降低了感知速率
  - 声波束按锥形方式传播,张开角约20～40度使得方向分辨率较差，聚焦性差，不能准确判断障碍物位置，方位角和形状精度差。
  - 软的物体表面将吸收大部分声音能量
  - 光滑的物体表面将形成镜面反射

##### 3.4.4.2 激光雷达

**原理**：

- 发射器用被校准的光束照亮目标
- 发射光束和接受光束同轴
- 机械装置控制旋转镜面实现2D或3D扫描

**优点**：
角度分辨率好，形状位置非常清晰，聚焦性好，可测角度范围大，可测距离远，测量的精度准确度高

**缺点**：
光速太快，不好测；激光穿透性差

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\3dd0b97f20d9dbd49ed79d3a3b0cda5.png" alt="3dd0b97f20d9dbd49ed79d3a3b0cda5" style="zoom:33%;" />

测量方法：

**时飞法（直接延迟时间测量法）**：和超声波一样
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\55a9071cc6b98b861170489a9db47f0.png" alt="55a9071cc6b98b861170489a9db47f0" style="zoom:50%;" />

**间接相位偏移测量法**：
发射器发射一个连续波。用具有不同频率的sin信号调制所携带信 号的波长。比较反射信号与所发送信号之间的相位差
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\ea589b9f4aa07ab78f3490fe852d3ec.png" alt="ea589b9f4aa07ab78f3490fe852d3ec" style="zoom: 67%;" />

##### 3.4.4.3 其他的：光学接近觉传感器

对置式

原理：用红外发光二极管发射红外光，根据接收器是否接收到红外光束确定是否有物体存在

回波式

原理：红外发光二极管发射红外光束，若物体接近到一定程度，红外光照射于物体上物体将反射红外光束，反射光束由接收透镜收集后汇集在光电二极管上，光电管有光电流输出

#### 3.4.5 其他传感器

##### 触觉传感器

##### 视觉传感器

- 原理：通过光学摄像机或红外、激光、超声、X射线对周围场景或物体进行探测成像
- 成像原理：小孔模型

## 4 机器人驱动

### 4.1 电机工作原理

注意：**力矩 X 转速 = 功率**

- 直流电机
  - 可以输出力矩和速度，如小车的直线运动、转弯等
  - 需要驱动芯片以及控制方式
- 舵机
  - 用于角度、位置伺服，如机械手转动
  - PWM波（占空比可变的方波）控制

### 4.2 无刷直流电机

- 无刷直流电机由**电机本体**、**位置传感器**、**电子换向电路**三大部分组成。 
- 电机主体由**主定子**、**主转子**组成。**主转子是永久磁铁**，**主定子是电枢**。当 定子绕组通直流电时，与转子作用产生电磁转矩，定子电流必须根据转 子的位置变化适时换向，才能获得单一方向的电磁转矩，使电机转动。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\0345198ce7849c9ab0dc7658b2bf01c.png" alt="0345198ce7849c9ab0dc7658b2bf01c" style="zoom:50%;" />

一般而言，无刷电机的绕组有星形联结方式和三角联结方式， 而三相星形联结（Y型）的二二导通方式最为常见

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\6f3aecdd45fd3c9079835d00f43648a.png" alt="6f3aecdd45fd3c9079835d00f43648a" style="zoom:50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\99ae30236a295ffa0da5b0e221c80ed.png" alt="99ae30236a295ffa0da5b0e221c80ed" style="zoom:50%;" />

三个重要的物理量：
电枢电动势Ea，电磁转矩T和电磁功率P；

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9c60c33c150e2bbef74531cc0b58e16.png" alt="9c60c33c150e2bbef74531cc0b58e16" style="zoom:50%;" />

转矩与转速的关系：
$$
U=E_a+I_a·R_a=K_en+I_aR_a\\
n=\frac{U-I_aR_a}{K_e}
$$
检查电机是否烧坏，可以通过测量电机绕组的阻值（短路会使电阻值变化）是否正常来判断。

#### 直流电机PWM调速

基本原理：动力源直接给电机供电；控制器输出的“小功率”控制信号；在控制信号作用下，电机驱动器把动力源的电压“调制”成“大功率”控制电压



 电路中串联**高速开关S**，把供给电机的连续电流离散化；电源供应的能量受到离散电流的调控；**当开关速度远大于电机响应速度时**，开关对电机转速平稳性的影响可以忽略不计，因此
设电机永远接通电源时，转速最大为Vmax，设占空比为D=t1/T，则有$$V_d=V_{max}·D$$

- 电机的转速与电机电枢电压成比例，而电机电枢电压与控制波形的占空比成正比； 
- 电机的速度与占空比成比例，占空比越大，电动机转得越快，当占空 比等于1，电机转速最大。

把连续变化的控制电压转化为固定频率的方波信号，方波占空比与电压大小成正比，利用方波信号控制电机 调速的技术称为脉宽调制(Pulse Width Modulation, PWM)控制技术。

PWM控制信号可以用MCU生成，通过程序调节占空比，无需使用模拟器件

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\b13071f6bf65d3190e834f3911cde71.png" alt="b13071f6bf65d3190e834f3911cde71" style="zoom:50%;" />

**带有正反转换向功能的调速原理**：

采用4个高速开关有序闭合、断开，可实现电机旋转换向

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\eb1db473a0e932653955b5db912ff12.png" alt="eb1db473a0e932653955b5db912ff12" style="zoom:50%;" />

### 4.3 有刷直流电机

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\8a35266f5d996813f189c62dab85fe4.png" alt="8a35266f5d996813f189c62dab85fe4" style="zoom:50%;" />

驱动芯片L289，双H桥式驱动器，可以方便的驱动两个直流电机。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\5eb3669a8dd7b3032d3f24d2da119d7.png" alt="5eb3669a8dd7b3032d3f24d2da119d7" style="zoom:50%;" />

#### 4.3.2 直流电机控制技术

**开环控制**：

- 利用电压与速度的对应关系，直接改变输入电压实现转速控制 
-  利用速度对时间的积分估计电机转过的角度，实现位置控制

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\91c13fa9084213e403580aa73a4cf63.png" alt="91c13fa9084213e403580aa73a4cf63" style="zoom:50%;" />

**闭环控制**：

 使用传感器检测速度（或位置）的实际输出值，以期望速度值与实际速度值之间的差值作为输入，经过控制器处理后生成电压控制信号，又称**闭环反馈控制或伺服控制**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\0df8900981a56c8815ab3dee3d22d4c.png" alt="0df8900981a56c8815ab3dee3d22d4c" style="zoom:50%;" />

优缺点：

- 优点：调速“刚性”好，精度高，响应快，抗扰动能力强 
- 缺点：系统复杂，可靠性低，成本高

在动态性能要求不高的场合，为降低系统应用的复杂度，可采用步进电机或舵机来替代伺服系统

#### 4.3.3 步进电机

- 步进电机通过脉冲信号进行控制，每输入一个脉冲信号，步进电机前进一步
- 输出转角与输入的脉冲的个数成线性关系，因此可以用开环控制的方式实现较为精准的 速度和位置控制

步距角： 电机每通过一个电脉冲，转子转过的角度θs 。常见的步距角为1.8°(0.9°,3.6°,7.5°)等

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\fe65458660e2668c91e2a8856bad5c8.png" alt="fe65458660e2668c91e2a8856bad5c8" style="zoom:50%;" />

**步进电机+步进驱动器+伺服控制器**优缺点：

优点：

- 可实现精密位置控制 
-  简单易用 
-  成本低

缺点：

- 存在“丢步”现象导 致精度下降
- 对负载波动适应性差

#### 4.3.4 减速比

传动比，也即减速比。指减速机构输入速度与输出速度之比，用“i”表示。即，i =输入速度/ 输出速度，并使输出力/力矩变为原来的i倍。

**用来减小速度，增大力矩**，功率 = 扭矩 × 转速

例：电机输入减速箱的速度1000n/min，输出速度10n/min，则减速比 i=1000/10=100 如电机输出力矩为Tin=0.1Nm，则输出力矩为Tout=Tin\*i=0.1Nm*100=10Nm

### 4.4 模拟舵机

减速齿轮组由电机驱动，其输出带到动一个线性的比例电 位器作为位置检测，该电位器把转角坐标转换为比例电压反 馈给控制线路板，控制线路板将其与输入的控制脉冲信号比 较，产生纠正脉冲，并驱动电机正/反转。

标准舵机有三条控制线，分别为电源线、地线和控制线。 控制线连接到控制芯片上。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c0a58d9d303ec4912d2e9303ce3cfdf.png" alt="c0a58d9d303ec4912d2e9303ce3cfdf" style="zoom:50%;" />

齿轮组：增大力矩，减小转速

控制电路：驱动器

比例电位器：检测角度

**舵机位置控制**

-  舵机转动角度由PWM（脉冲宽度调制）信号的占空 比来实现；
-  PWM周期为20ms，脉宽分布在0.5～2.5ms之间
- 不同脉宽对应不同转角位置。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c4933ed510c4f0880a20e9616ef786b.png" alt="c4933ed510c4f0880a20e9616ef786b" style="zoom:50%;" />

##### 电机减速器

- 传递能量传递：驱动----->机构
- 改变运动速度
- 改变运动力量
- 改变运动方向

### 4.5 气动设计

一个简易的气动系统由：气压元件、控制元件、执行元件和 辅助元件组成。

**方向控制回路**

**单作用气缸换向回路**

掌握几位几通的概念：

- 几位：看有几种回路可以切换
- 几通：看每一个回路里面有几个连接处

位就是有几个工作状态，几个框就是几位。

通就是与边框交点，取一个小框，线与边框交点是多少就是几通。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\20a8d3f98852bfdec769681e689ad62.png" alt="20a8d3f98852bfdec769681e689ad62" style="zoom:67%;" />

对于a，两个框，一个框3个交点，两位三通阀。

**双作用气缸换向回路**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\91e47599eb864b850e547de57a539f3.png" alt="91e47599eb864b850e547de57a539f3" style="zoom:67%;" />

在电磁阀作用前后，电磁阀控制的气缸，从一个位转换成另一个位，其左右相连的线路顺延继承给下一个。以c为例子，上边按钮没按之前，向有杆腔输气，杆左移。按下后，对上面气缸，换位，左侧变为从右箭头处输入，右侧变为从右箭头处输出，输出过来的气体使得下方的二位五通气缸换位，这样变为从右箭头处输入气体并输出到无杆腔处，杆右移。

**气动肌肉**

气动人工肌肉是一种体积小巧、柔软、重量轻、工作简单、容易控制的仿生学产品。它只能收缩到一定长度，具有较好的安全性

**液压驱动**

- **将液体压力转化为机械能**
- 利用**不可压缩的流体**，将作用于某一点的力传递到另一点，这种流体通常是**工业液压油** 
- 由液压源、伺服阀、传感器、执行机构等构成
  - 优点： 重量轻、尺寸小、动作平稳、快速性好、产生的力 力矩非常大 。
  - 缺点： 易漏油、维护困难；不确定性和非线性因素多，控制和校正不如电气式方便

**新型驱动**

气体软体机器人

### 4.6 齿轮

#### 齿轮基本参数：

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\be134e49618cb9537ec3d17f5a1704d.png" alt="be134e49618cb9537ec3d17f5a1704d" style="zoom: 50%;" />

模数m：人为地把$$\frac{p_i}{\pi}$$规定为一些简单的有理数，该比值称为模数，用m表示。模数越大，齿厚就越大，齿轮的承载能力就越高。

分度圆d：是齿轮上一个人为地约定的轮齿计算的基准圆，规定分度圆上的模数和压力角为标准值。分度圆又称节圆 。国际压力角标准值20°。

分度圆直径$$d=mz$$
基圆直径$$d_b=dcos\alpha=mzcos\alpha$$
基圆齿距$$P_b=\pi d_b/z=\pi mcos\alpha$$
齿顶圆直径$$d_a=d+2h_a$$
齿根圆直径$$d_f=d-2h_f$$
法结$$P_n=P_b=\pi mcos\alpha$$

欲使两齿轮正确啮合，两轮的模数必须相等

#### 定常传动比

对齿轮传动的基本要求是$$i_{12}=\omega_1/\omega_2=C（常数）$$

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\af2e8239a8a94f69d414ed44c3f92cc.png" alt="af2e8239a8a94f69d414ed44c3f92cc" style="zoom: 50%;" />

#### 减速比

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\651896f6bc1aa4aea99ff33429403ae.png" alt="651896f6bc1aa4aea99ff33429403ae" style="zoom:50%;" />

定轴轮系传动比$$i_{首末}=\frac{\omega_首}{\omega_末}=\frac{z_末}{z_首}$$

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\06e374b011adb55fc8531c292c23622.png" alt="06e374b011adb55fc8531c292c23622" style="zoom: 67%;" />

转向关系，外啮合反向，内啮合同向

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c3116042002266a59326f0feee025c3.png" alt="c3116042002266a59326f0feee025c3" style="zoom: 50%;" />

注意有可能一个轴连两个z的齿轮，计算传动比时要注意是哪个齿轮啮合的。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\291543fc38edfe900074b1e85cc9940.png" alt="291543fc38edfe900074b1e85cc9940" style="zoom:50%;" />

### 4.7 连杆

铰链四连杆是平面四杆机构的基本形式，其他形式可以认为是它的演化形式。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\340e20e6c0905b2ae633fe931fd60ef.png" alt="340e20e6c0905b2ae633fe931fd60ef" style="zoom:50%;" />

原动件AB的运动经过一个不直接与机架相联的中间构件BC才能传动从动件CD

**曲柄摇杆机构**

曲柄摇杆机构：铰链四杆机构的两个连架杆中， 有一个为曲柄，另一个为摇杆。

曲柄摇杆机构的条件：平面四杆机构具有整转副👉则可能存在曲柄

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\02e85e08e37ecc1bbf08c61b654acd3.png" alt="02e85e08e37ecc1bbf08c61b654acd3" style="zoom:50%;" />

**双曲柄机构**

 若铰链四杆机构中的两个连杆架均为曲柄，则称其为双曲柄机构。

**双摇杆机构**

若铰链四杆机构中的两个连杆架均为摇杆，则称其为 双摇杆机构。

**导杆机构**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\81efc9b5f01a34e8c85ff0a179636c7.png" alt="81efc9b5f01a34e8c85ff0a179636c7" style="zoom:50%;" />

**曲柄滑块机构**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\ea4070d4def9782665cf984f25f4046.png" alt="ea4070d4def9782665cf984f25f4046" style="zoom:50%;" />

**曲柄滑块结构的位移与速度分析**
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7e910daf9a10884396b4354de6782c0.png" alt="7e910daf9a10884396b4354de6782c0" style="zoom:50%;" />已知若干长度、角度、角速度

**位移分析**：用四边形的矢量方程，实部与虚部求解

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\be0eda5b87493ef72f89093637a61f8.png" alt="be0eda5b87493ef72f89093637a61f8" style="zoom:50%;" />

**速度分析**：对位移矢量方程求时间导，再实部与虚部相等。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\bd54ade637a565a22dfd2e915ec86fc.png" alt="bd54ade637a565a22dfd2e915ec86fc" style="zoom:50%;" />

**平行四边形机构**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\13f3a8997270cdc60cef00545beb021.png" alt="13f3a8997270cdc60cef00545beb021" style="zoom:50%;" />

### 4.8 轴承

按照承受的载荷方向不同，滚动轴承分三种

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9f9cab9bf6bd33c595f9d86f3bbd02b.png" alt="9f9cab9bf6bd33c595f9d86f3bbd02b" style="zoom:50%;" />

接触角a：滚动体的荷载方向线与轴承径向平面之间的夹角。a越大，可以承受的轴向力越大。

**轴承的安装**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\4663492f48646aeca04655fce7cafe5.png" alt="4663492f48646aeca04655fce7cafe5" style="zoom:50%;" />

### 4.9 机器人运动学

运动学:是指机器人连杆的位置和姿态（简称：位姿）与关节角度关系的理论

正运动学：已知关节角，求连杆末端的位姿
逆运动学：已知连杆末端的位姿，求关节角度

正运动学一般用于检验逆运动学是否正确，逆运动学难度高于正运动学。

运动学方程是非线性的，很难得到逆运动学的解，可能无解或者多解。

**转动特性**：

最基本的转动是分别按照x、y、z轴转动，分别称作滚动，俯仰，偏摆

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\b50c23456460c3ce8b3c9adcc3ea09d.png" alt="b50c23456460c3ce8b3c9adcc3ea09d" style="zoom:50%;" />

对应的**旋转矩阵**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\4f81713703492c1d4451d5b3c204222.png" alt="4f81713703492c1d4451d5b3c204222" style="zoom:50%;" />

**转动特性**，三个矩阵的逆序依次乘法，按照zyx的顺序

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\577a59a0baaa3b3631f1ee2bda63b6b.png" alt="577a59a0baaa3b3631f1ee2bda63b6b" style="zoom:50%;" />

## 5 机器视觉

### 5.1

被动视觉传感器：被动接受光，如摄像头

主动视觉传感器：主动发射光，如深度相机

一对被动视觉传感器构双目相机，获得深度，公式为：
$$
d=\frac{c}{2tan(a/2)}
$$
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\ec8a7eeb0dbc7b508b8a5f320e4f8e8.png" alt="ec8a7eeb0dbc7b508b8a5f320e4f8e8" style="zoom:50%;" />

距离越远误差越大。

### 5.2 应用场景

瑕疵检测、土壤分析、表计读数、文字识别、病灶检测、人脸识别、手 势分析、疲劳驾驶、工件建模、地形测绘、视觉导航、视觉定位、车辆 预测、场景解析、虚拟渲染、头手跟踪、娱乐等

视觉提供了一种**几何测量**的工具，也提供了一种**语义认知**的工具。各种视觉工具是两种工具的组合。

### 5.3 机器人视觉系统

从自然界到信息的通路在机器人上的实现

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\ca081e426ca675ba495094fda49950c.png" alt="ca081e426ca675ba495094fda49950c" style="zoom:50%;" />

### 5.4 图象建模

图像是定义在CCD阵列下的离散函数。
$$
I:(u,v)\in[0,W-1]\times[0,H-1]→q\in R^N\\
q=I(x)
$$
做个图像分类的都干过，一个图像分为若干个像素点，每个像素点有RGB三原色，存放着一个颜色元素。

### 5.5 相机建模

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\54ce5bf45afe9af6bdda03012dbe07e.png" alt="54ce5bf45afe9af6bdda03012dbe07e" style="zoom:50%;" />

**内参**

焦距（fx, fy）

- 表示光学中心到图像传感器的距离，以像素为单位。
- fx 和 fy 分别表示图像在水平和竖直方向的放大倍数（像素/mm）。
- 一般情况下 fx ≈ fy。

主点（光学中心）（cx, cy）

- 图像坐标系的原点位置，通常位于图像中心附近。

像素倾斜系数（s）（在矩阵第一行第二列处）

- 通常为0，表示传感器的像素横竖完全垂直。

内参只与相机本身有关

**外参**
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\33f4a5951a098a36d2f59c28f994479.png" alt="33f4a5951a098a36d2f59c28f994479" style="zoom: 80%;" />

其中R是旋转矩阵，表示相机相对于世界坐标系的旋转关系（姿态）。
t是平移向量，表示相机坐标系原点相对于世界坐标系原点的位置。

**变量**

最右侧的P其实是物体在世界坐标系的坐标列向量，最下边加个1是为了与平移向量吻合。

**前向计算过程**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\24c31c2cdcf0af73cd217fe44519273.png" alt="24c31c2cdcf0af73cd217fe44519273" style="zoom:67%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\734183068aec23e3fe41420df21a673.png" alt="734183068aec23e3fe41420df21a673" style="zoom:67%;" />

#### 镜头畸变

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\a575496cb7bc40fe9b011d653a7ce8f.png" alt="a575496cb7bc40fe9b011d653a7ce8f" style="zoom:50%;" />

**内参估算**

相机标定，借用外部已知尺寸的物体，解算出内参。

张正友标定法：

- 采用棋盘格作为已知尺寸的物体，利用平面特性方便求解  
- 棋盘格的角点检测相对简单，可靠性高  

**外参应用**

- **AR的应用原理**就是外参的应用

- 基于指定尺寸平面，可以估计出平面和相机的外参，也就是相机在世界坐标系下的位姿
- 如果在世界坐标系下，增加一个虚拟点，可以计算出在图像中的成像  

## 6 机器人定位

VSLAM 即 Visual Simultaneous Localization and Mapping，主要是指如何用相机解决定位和建图问题。当用相机作为传感器时，我们要做的，就是根据一张张连续运动的图像(它们形成一段视频)，从中推断相机的运动，以及周围环境的情况。

**SLAM** (`simultaneous localization and mapping),也称为CML (Concurrent Mapping and Localization`), 即时定位与地图构建，或并发建图与定位。

定位问题：确定机器人在世界（全局）坐标系中的位置/位姿  

全局坐标系：世界坐标系
局部坐标系：相对于世界坐标系某一部分的相对坐标系

### 6.1 里程估计

- 根据传感器感知信息推导机器人位姿（位置和角度）变化  
- 用途：
  -  航位推算 (Dead-reckoning)：基于已知位置， 利用里程估计，推算现在位置

里程估计方法  

- 基于机器人运动感知信息，结合运动学模型
  - 电机码盘
  - IMU（惯性单元，加速度计+陀螺仪）
- 基于环境感知传感器信息，通过最佳匹配估计
  -  激光里程计
  -  视觉里程计（VO）

#### 6.1.1 轮式里程估计

 轮式里程计是典型的里程估计方法，但含有多种误差



<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\080cbe10cb776aa35595b0ef9d5ae6d.png" alt="080cbe10cb776aa35595b0ef9d5ae6d" style="zoom:50%;" />

此为基于电子码盘的轮式移动机器人里程估计。
（1）根据电子码盘获得轮子转速，对应第四个公式，n是码盘测得的电机转速（转/分），分母是齿轮减速比

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\bc8955d2dcc63828f53f9ed0a3fdc0d.png" alt="bc8955d2dcc63828f53f9ed0a3fdc0d" style="zoom:50%;" />

（2）结合运动学模型计算参考点速度

（3）假设短时间内为匀速运动，计算位姿变化（误差会积累）

系统误差：

- 轮半径
- 轮子安装精度误差（不平行，两边距离不相等）
- 编码器精度误差
- 采样精度误差
- 齿轮减速比精度

偶然误差：

- 地面不平
- 轮子打滑

**基于姿态变化的航位推算**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\00a99beabe82217a51ac8e29915d57c.png" alt="00a99beabe82217a51ac8e29915d57c" style="zoom:50%;" />

**基于惯性单元的里程估计**  

惯性单元IMU  

- 一般含有三轴的加速度计和三轴的陀螺仪
- yaw不准

- 优点
  - 全天候
  - 采样频率高
  - 短时精度较好
- 缺点
  - 随着时间的增长累积误差较大，无法满足移动机器人长距离精确定位的要求，需要融合其它传感器进行组合导航

#### 6.1.2 激光里程计

基本原理：

采用ICP(Iterative Closest Point)算法  

- 估计P’集合点与P集合点的初始位姿关系  
  - P'是P点相对于里程计的坐标系下位置改变后的坐标
- 根据最近邻域规则建立P’集合点与P集合点的关联  
  - 构建旋转矩阵，使得P'旋转之后成为P点
- 利用线性代数/非线性优化的方式估计旋转平移量  
  - 优化旋转矩阵
- 对点集合P’的点进行旋转平移  
- 如果旋转平移后重新关联的均方差小于阈值， 结束
- 否则迭代重复上述步骤  

在开始点，用激光构建周围障碍的点集合；在终止点，也用激光构建周围障碍点的集合。将这终止点集合平移旋转直到两个集合大致吻合，集合的相对位姿变化便是两个点的相对位姿变化。

### 6.2 定位问题

#### 6.2.1 GPS

GPS基本原理：多颗卫星的纯距离测量的融合

由**空间端、 控制端和用户端**三部分组成，也称为GNSS

- 空间端由卫星组成
- 控制端由地面天线站等构成，主要作用监测控制卫星的运行，并对卫星进行时间同步  
- 用户端，采用三点法进行定位

误差：（多径时需要更多卫星保证可靠性）

遮挡问题：室内难以使用

多路径问题：楼房太高导致最短路径变成反射后的一条路径

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\2e45c333f96d287c20ff9812157800d.png" alt="2e45c333f96d287c20ff9812157800d" style="zoom:50%;" />

#### 6.2.2 全局视觉观测定位

利用外部视觉模拟GPS，对机器人进行测量定位

搭建一套外部视觉系统，识别机器人，确定其位置  

与GPS相同，在内部搭建多个相机实现定位

应用约束：

- 摄像头有视野范围约束，当环境较大时需要多个摄像头
- 机器人上也需要有一定的标识方便图像识别定位

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\057fbc7951788ffb4bf72f03d341a70.png" alt="057fbc7951788ffb4bf72f03d341a70" style="zoom:50%;" />

#### 6.2.3 通过本体视觉识别空间标识

**基于环境人工标识的定位**

在环境中部署特殊标签，降低成本，确保可靠性，如

- 地面二维码：适合仓库运输机器人
- 高反射物质：适合叉车

**基于环境自然标识的定位**  

路标

**基于空间标识的定位原理**  

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9f0aa1215c13be8a31dbdd9a85af04a.png" alt="9f0aa1215c13be8a31dbdd9a85af04a" style="zoom:50%;" />

已知三个星星的坐标，并通过传感器测得r或者θ，进行求解坐标

#### 6.2.4 概率融合全局与里程观测

定义：在给定环境地图(对环境的已有知识)的条件下，根据机器人的运动控制/里程估计信息和传感器感知数据(当前认知)估计机器人相对于环境地图的坐标

能够避免非唯一数据关联引起的歧义

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\c862e3ee625f87282ac22cdb0a8cdad.png" alt="c862e3ee625f87282ac22cdb0a8cdad" style="zoom:50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\9d18aab378a8c1923f7a0492e02b684.png" alt="9d18aab378a8c1923f7a0492e02b684" style="zoom: 67%;" />

例子：机器人在一扇门前，但是不知道在哪扇门。

黑色是对由里程计对自己位置的估计，起始没有传感器信号，机器人认为是均匀分布。红色是传感器得到的信息，得到三个门的信息，在三个门处得到了三个峰值。把这两个概率密度函数融合，融合完后是机器人经过一次传感器信息后对自己位置的判断，还是三个峰值。然后机器人向右从第一扇门移动到第二扇门，根据里程计得到的信号，机器人对自己的位置判断在概率分布上整体向右平移（第一条），此处传感器再次得到三个峰值的信息（第二条）。然后再将这两个概率融合，是机器人再经过一次传感器信号后对自己位置的判断，此时在第二扇门处出现了峰值（第三条），说明机器人在第二扇门处。然后机器人再向右移，再根据里程计得到相同的右移（第四条），唯一大峰值也右移，这样即使机器人不在门前也可以知道自己的位置了。

## 7 机器人规划

### 7.1 轨迹生成 

**生成光滑一维轨迹**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\6d6fb220aee07aa2285e1be83292cac.png" alt="6d6fb220aee07aa2285e1be83292cac" style="zoom:50%;" />

五次多项式轨迹$$x(t)=p_5x^5+p_4x^4+p_3x^3+p_2x^2+p_1x+p_0$$

约束条件：<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\195db71975edb8f637f6760d5881c28.png" alt="195db71975edb8f637f6760d5881c28" style="zoom:67%;" />

列出矩阵求解：

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\2aaa748a0cbfceb7d926443d1c186ff.png" alt="2aaa748a0cbfceb7d926443d1c186ff" style="zoom:67%;" />

**生成光滑多段轨迹**

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\caffafb2f1427d547e113a6bf4883c4.png" alt="caffafb2f1427d547e113a6bf4883c4" style="zoom: 50%;" />

五次多项式轨迹$$x(t)=p_5x^5+p_4x^4+p_3x^3+p_2x^2+p_1x+p_0$$

约束条件：<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\8f86d28957254bb35f146748c0cb2aa.png" alt="8f86d28957254bb35f146748c0cb2aa" style="zoom: 67%;" />

列出矩阵求解：

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\8998ba0e1befa85a8758538533d7802.png" alt="8998ba0e1befa85a8758538533d7802" style="zoom: 80%;" />

### 7.2 路径搜索

基于采样的方法

- PRM(Probabilistic roadmap)
- RRT(Rapidly exploring random tree)
- RRT*

基于搜索的方法

- 图搜索：DF，BFS
- dijkstra，A*
- Jump Point Search

#### PRM

基于概率采样的路径

- 均匀生成采样点
- 将与障碍物接触的点给删除
- 领域点计算：在距离为r的圆内均为领域点，将其连接
- 碰撞检测：连线是否与障碍物相交

优点：产生的roadmap可以被复用

缺点：对于给定的起点和终点，非最短路径，效率低

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\189ef107d6a5fa778807cda3706ac74.png" alt="189ef107d6a5fa778807cda3706ac74" style="zoom:50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\4059a14fe1a90cdfa54cb07a04a45d8.png" alt="4059a14fe1a90cdfa54cb07a04a45d8" style="zoom: 67%;" />

#### RRT

从根或者已生成树开始，规定两点之间的步长固定（也可不固定）

- 在地图上随机采样点
- 采样点与最近点相连接，生成树
- 如果不经过障碍物，此路径满足条件，可以加上

优点：容易添加对目标点的引导，效率增加

缺点：无法删除已生成树，但不一定是最短，没有纠错机制

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\f2658a953e7b936d58eef3ae4df08b2.png" alt="f2658a953e7b936d58eef3ae4df08b2" style="zoom:67%;" />

#### RRT*

在添加新路径时，面对若干满足的采样点，RRT*会判断从起点到哪个采样点路径总距离最短。

会改变过去的信息，增加了Rewrie（重连）函数，即在采样之后与最短路径连接后，考虑在某一个定长的圆的范围内，其内的点是否可以连接到新采样的点（用到初始点的距离进行判断），可能改变一些点的父节点。

#### Dijkstra

加了权重的广度优先算法

- **优点**: 
  - **准确性**: 总是能找到最短路径。
  - **简单性**: 实现相对简单。
- **缺点**: 
  - **效率较低**: 算法需要遍历图中的大多数节点，可能导致较高的计算成本。
  - **实时性差**: 在动态环境中可能不适用，因为它不能快速适应环境的变化

#### A*

在Dijkstra's算法的基础上，加了对于距离目标点的预测方向，因而有了更强的目的性（启发式算法）

- **优点**: 
  - **效率较高**: 通过启发式函数，A*算法能够减少需要遍历的节点数量，从而显著提高搜索效率。
  - **准确性**: 同样能够找到最短路径。
  - **灵活性**: 通过改变启发式函数，可以很容易地对算法进行定制。
- **缺点**: 
  - **启发式函数选择**: 启发式函数的选择对算法的性能有很大影响，不恰当的启发式函数可能导致搜索效率降低。
  - **实时性**: 虽然比Dijkstra算法更高效，但在高度动态的环境中仍可能面临挑战。

#### JPS

跳点寻路算法，JPS是在A*算法上的一种改进，算法在遍历邻居的时候不是直接将该点加入到队列中，而是将该方向的跳点加入到队列中（关于跳点的判断是整个JPS算法的核心）。这样整个搜索队列仅保存少量的点，可以以非常好的优化搜索时间。

三者都是最优解

## 8 机器人集群

### 8.1 基于Virtual Structures的编队控制

核心思想：

- 集群编队结构表示：将整个集群用virtual structures表示为一个世界坐标系下的整体（virtual rigid body）（VRB）
- 多目标需求：基于势场法表示集群中每架无人机编队保持、相互躲避、障碍物避障的需求。
- 控制：在VRB坐标系下统一上述各个势场得到相应的控制指令。

基于VRB的集群编队控制

VRB坐标：平移+旋转$$\delta_v(t)=(p_v(t),R_v(t))$$
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\67937f32450604df7cdfd1edbc0702a.png" alt="67937f32450604df7cdfd1edbc0702a" style="zoom:67%;" />
<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\7e7ea9c1b22287157c1bc070fd5883d.png" alt="7e7ea9c1b22287157c1bc070fd5883d" style="zoom:67%;" />

期望状态就是平移加旋转。

对时变编队，r(t)是t的函数，求导也要连带求导。

对固定编队，r(t)是常数，求导时视为常数。

####  在障碍物环境下基于势场法的编队控制

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\66c702c93a867ba6d38b58fe480c425.png" alt="66c702c93a867ba6d38b58fe480c425" style="zoom: 50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\771f1a7b87f842b70baa67e51275e44.png" alt="771f1a7b87f842b70baa67e51275e44" style="zoom: 50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\607b9220854268c6b3980ee0a181ac2.png" alt="607b9220854268c6b3980ee0a181ac2" style="zoom:50%;" />

### 8.2 基于速度障碍物的多智能体避障算法

VO（velocity obstacle）

问题简化：直径2r的机器人（红色）避开一个障碍物（蓝色）等价于把机器人看作质点，障碍物膨胀r。 VO的直观概念：对于红色机器人，任何落在浅蓝色区域内速度矢量被称作速度障碍物（VO），因为这 些速度将最终导致与障碍物碰撞。

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\30e1b7df7ea6236488f3602a7c7d671.png" alt="30e1b7df7ea6236488f3602a7c7d671" style="zoom:50%;" />

针对移动障碍物的VO. 

VO 𝐴|𝐵(𝒗𝑩)定义为：A将与速度为𝒗𝑩的B相撞的速度空间

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\a8690d060f49c3ea90b9c0b3ee96578.png" alt="a8690d060f49c3ea90b9c0b3ee96578" style="zoom:50%;" />

针对多个移动障碍物的VO

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\5c6820394cb923e49d905385f93af28.png" alt="5c6820394cb923e49d905385f93af28" style="zoom:50%;" />

<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\5a815185306efe7a2c7a95737e5c388.png" alt="5a815185306efe7a2c7a95737e5c388" style="zoom:50%;" />

缺点：VO的振荡问题，根本原因：：每个机器人只考虑其它机器人当前的速度， 而不考虑其他机器人下一个控制周期的速度。

优点：复杂度低

#### 改进

RVO：reciprocal velocity obstacle

$$RVO_{A|B}(v_A,v_B)$$定义为速度VA的A将于速度为VB的B相撞的速度空间

基本思想：假设所有个体都运行RVO，选择速度$$v_A'$$，它是个体A的当前速度$$v_A$$和VO区域外的一个安全速度v的平均，即所有个体都适当减小一些速度的改变量$$||v_A'-v_A||$$
