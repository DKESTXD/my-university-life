<img src="D:\新建文件夹 (3)\WeChat Files\wxid_uzqin1xmai4h22\FileStorage\Temp\676d04f611e452a287b118125bb151c.png" alt="676d04f611e452a287b118125bb151c" style="zoom:67%;" />

人声与伴奏除了音色上的显著差别之外，在频率范围上也有明显差别，

如人声为80-1.5KHZ，乐器为20-6KHZ，而且不同的音乐人声不一样，乐器也不一样，因此能否抓住其频率的差异，用滤波器提取特定频率范围的音频内容。

因为伴奏是从始至终一直有的，而人声是相当于嵌入进去的，  所以可以先专门提取人声，再在原音频中剔除人声得到伴奏（实际上这个方法不现实，原因后面再说）。

## 频率分离法

#### 只用滤波器

试一试：用一个带通滤波器取3000-20000HZ的范围为人声，取小于200与大于20000为伴奏

```matlab
% 加载音频
[x,fs]=audioread('./original audio/someone.wav');

%% 1. 人声提取（带通滤波器）
bpFilt=designfilt('bandpassfir','FilterOrder',200,'CutoffFrequency1',3000,'CutoffFrequency2',20000,'SampleRate', fs);
x_vocal=filter(bpFilt, x);
x_vocal=x_vocal*8;
%% 2. 伴奏提取（低通 + 高通滤波器组合）
% 低通滤波器：鼓/贝斯等低频成分
lpFilt=designfilt('lowpassfir','FilterOrder',200,'CutoffFrequency',200,'SampleRate',fs);
x_low=filter(lpFilt, x);

% 高通滤波器：合成器/弦乐高频
hpFilt=designfilt('highpassfir','FilterOrder',200,'CutoffFrequency',20000,'SampleRate',fs);
x_high=filter(hpFilt, x);

% 合成伴奏部分（非人声）
x_music=x_low+x_high;

%% 3. 保存结果
audiowrite('vocal_only.wav',x_vocal, fs);
audiowrite('music_only.wav',x_music, fs);
Omegaplot(x,fs);
Omegaplot(x_vocal,fs);
Omegaplot(x_music,fs);
```

效果不用看 肯定相当之差

对hundun来说，人声频率范围取2000-20000HZ，伴奏为<400，>15000HZ

这是因为人声与伴奏肯定会有重叠的地方，即使真的存在只有人声的频率范围，其提取出来的人声也不完整，且在剩余频率范围内存在较大残留，如果再修改频率范围来寻找只有伴奏而无人声的频率段，提取出的伴奏也不完整。因此人声与伴奏频率重叠部分的分离是难点。

再者，上述频率范围的确定，是手工一次次试出来的，而且针对不同的歌曲，对应的频率范围也不一样，这种一刀切的方法不光在决定在何处切的时候麻烦，而且切出的结果也不完整

#### 第一次升级：加掩码

**频率掩码**：即指定一段频域上对应的值为0或1，相当于掩盖住了一部分频域，这种称之为硬掩码，本质上还是滤波器。与之对应的软掩码（高斯掩码），可以将一段频域上的值乘一定的比例，使之连续，比较平滑。

```matlab
clear;close all;clc;
%% 读取音频
[x,fs]=audioread('./original audio/someone.wav');
x=mean(x,2);  % 单声道处理

%% Step 1: 初步滤波（双通带滤波器）
bpFilt=designfilt('bandpassiir','FilterOrder',200,...
    'HalfPowerFrequency1',3000,'HalfPowerFrequency2',20000,...
    'SampleRate',fs);
% 滤波器提取人声频段
vocal_pre=filter(bpFilt,x);

%% Step 2: STFT 分析（用于动态掩码增强）
nfft=1024;
hop=256;
win=hamming(nfft,'periodic');
[S,f,t]=stft(x,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);
powerSpec=abs(S).^2;

% 归一化功率谱
powerNorm=powerSpec./(max(powerSpec(:))+eps);

% 构造频域高斯掩码
f0=10000;
sigma=9000;
gaussian_mask=exp(-((f-f0).^2)/(2*sigma^2));
gaussian_mask=repmat(gaussian_mask,1,length(t));

% 动态谱掩码
vocal_mask=gaussian_mask.*powerNorm;

% 平滑掩码
vocal_mask=imgaussfilt(vocal_mask,1);

%% Step 3: ISTFT 重建
S_full=stft(vocal_pre,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);
S_vocal=vocal_mask.*S_full;
S_music=S_full-S_vocal;

x_vocal=istft(S_vocal,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);
x_music=istft(S_music,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);

x_vocal=real(x_vocal);
x_music=real(x_music);
%% Step 4: 保存音频
audiowrite('vocal_filtered_masked.wav',x_vocal,fs);
audiowrite('music_filtered_masked.wav',x_music,fs);

%% Step 5: 可视化掩码与频谱
figure;
subplot(3,1,1);
imagesc(t,f,20*log10(abs(S)+eps));axis xy;
title('原始频谱');xlabel('时间 / s');ylabel('频率 / Hz');colorbar;

subplot(3,1,2);
imagesc(t,f,vocal_mask);axis xy;
title('频谱增强掩码');xlabel('时间 / s');ylabel('频率 / Hz');colorbar;

subplot(3,1,3);
imagesc(t,f,20*log10(abs(S_vocal)+eps));axis xy;
title('增强后人声频谱');xlabel('时间 / s');ylabel('频率 / Hz');colorbar;
```

STFT短时傅里叶变换分析

ISTFT逆变换还原

谱熵法估计频率中心（失败，测得的中心频率总是偏低）





提取出的人声相比于上一个，确实更丝滑了（这里出了个bug，研究半天没明白，ai就算了，他连频率都分不清，跑完代码本应生成人声的文件没声音，生成的伴奏确是效果比较好的人声？）；

#### 第二次升级：相位差谱

一般的音乐音频，都是有两个通道的，左声道和右声道，然而左声道与右声道并不是完全相同的， 在录歌的过程中，一般人是正对麦克风的，左右通道接收到的人声是同相位的，而对于伴奏，有的是在录音室中直接播放的，这样就会导致伴奏的声音会从不同方向传入麦克风，导致左右通道产生相位差。但是现在很多都是直接数字混音加进去的，这就导致人声和伴奏都是同相位的，也分不开。甚至还有人声是不同相位的，伴奏是同相位的情况。

可以构建两个通道的相位差谱，再用高斯掩码，保留相位差为零的部分，从某种程度上可以让其自行确定频率范围

```matlab
%% 读取立体声音频
[x,fs]=audioread('./original audio/hundun.wav');  % 假设是双声道
x1=x(:,1);  % 左通道
x2=x(:,2);  % 右通道

%% STFT 设置
nfft=1024;
hop=256;
win=hamming(nfft,'periodic');

[S1,f,t]=stft(x1,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);
[S2,~,~]=stft(x2,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);

%% 相位差谱
phase_diff=angle(S2.*conj(S1));  % ∠(S2 × S1*) ∈ [-π, π]

% 计算相位差稳定性（小波动 = 人声方向）
mask=exp(-(phase_diff.^2)/(2*0.1^2));  % 高斯掩码，σ控制方向性宽度

%% 应用掩码增强人声
S1_masked=S1.*mask;

%% ISTFT 重建
x_vocal=istft(S1_masked,fs,'Window',win,'OverlapLength',nfft-hop,'FFTLength',nfft);
x_vocal=real(x_vocal);

%% 保存结果
audiowrite('vocal_phasediff.wav',x_vocal,fs);

%% 可视化
figure;
subplot(2,1,1);
imagesc(t,f,abs(S1)); axis xy; title('原始 STFT（左通道）'); xlabel('时间 / s'); ylabel('频率 / Hz');

subplot(2,1,2);
imagesc(t,f,mask); axis xy; title('相位差掩码'); xlabel('时间 / s'); ylabel('频率 / Hz');
```

对于someone这段音频，提取人声效果还行，但是人声音色存在失真

于是将掩码平滑处理一下

对于hundun这段音频，这就属于是伴奏同相位，但是人声不同相位的。分离出来的部分是伴奏，且较为不完整。

这个方法还是有很大的局限性的（受限于人声和伴奏的混音方式），但是效果确实比之前好一些了，尤其是其可以自己决定滤除的频率。



说实话，这个办法的核心还是频率分离的方法，只是把前者传统的直接截取一大段频率范围的办法，转变为截取若干个小段，并且截取的依据是左右通道的相位差。





从这我们看出，只从频域的角度来看，是无法分开那些同频率的人声和伴奏的，因为他们在频域上已经只代表一个幅值，已经是一个整体了，matlab只能分析不同频率的能量高低，却不能分析是谁的能量，即不能分别音色，这是核心问题，也是致命问题。

## 中通道法

继承上述左右声道的思想，立体声混音还有一个特点，就是人声一般居中，在通道中间，即人声相当于左右声道的共有部分（相位差为0）；而伴奏是左右声道的差异部分（存在相位差）。根据这个想法，我们也可以实现一下加强人声与削弱人声。（实际上也与人声混音有关）

这种也是卡拉OK常用提取伴奏的方法

```matlab
[x,fs]=audioread('./original audio/rolling.wav'); % x: N × 2 stereo
L=x(:,1);
R=x(:,2);
% 中通道
center=(L+R)/2;
% 边通道
side=(L-R)/2;
audiowrite('vocal_only.wav',center/max(abs(center)),fs);
audiowrite('accompaniment_only.wav',side/max(abs(side)),fs);

```

对于这种方法，边通道的效果比中通道更好，因为边通道求两者差值，直接减去了相同的部分，而中通道只是加强了中间部分，而并没有滤除差异部分，所以听感上只是人声加强了。

这种其实也和音乐本身有较大关系，在某些音乐上有奇效，但是在大部分音乐上只是起到一个削弱人声的作用，并没有完全消除。

someone的伴奏提取效果一般，但是能听出来人声的削弱

rolling伴奏提取效果非常好

hundun的人声相位情况和伴奏是反过来的，因此提取差值提取的是人声，效果也非常好。









