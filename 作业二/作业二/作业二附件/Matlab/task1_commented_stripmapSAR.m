%%================================================================
%% 文件名: task1_commented_stripmapSAR.m
%% 项目: 条带式SAR仿真 —— 点目标回波生成与距离多普勒(RD)算法成像
%% 作业一：逐行注释原始程序
%%================================================================
clear;clc;close all;  % 清除工作区变量、清空命令窗口、关闭所有图形窗口

%%================================================================
%% 常量参数
C = 3e8;                    % 电磁波传播速度（光速），单位 m/s

%% 雷达系统参数
Fc = 1e9;                   % 载波频率，1 GHz
lambda = C/Fc;              % 载波波长 λ = c/f，单位 m

%% 目标区域参数
Xmin = 0;                   % 成像区域方位向最小坐标（m）
Xmax = 50;                  % 成像区域方位向最大坐标（m）
Yc = 10000;                 % 成像区域中心的地距（m）
Y0 = 500;                   % 成像区域距离向半宽度（m），成像宽度为 2*Y0

%% 轨道参数
V = 100;                    % SAR平台飞行速度，单位 m/s
H = 5000;                   % 平台飞行高度，单位 m
R0 = sqrt(Yc^2 + H^2);     % 平台到成像中心的最短斜距（m）

%% 天线参数
D = 4;                      % 方位向天线物理长度（m）
Lsar = lambda*R0/D;         % 合成孔径长度（m），由分辨率要求决定
Tsar = Lsar/V;              % 合成孔径照射时间（s），即目标在波束内的时间

%%================================================================
%% 慢时间域（方位向）参数设计
Ka = -2*V^2/lambda/R0;     % 多普勒调频率（Hz/s），负号表示频率随时间减小
Ba = abs(Ka*Tsar);          % 多普勒信号带宽（Hz）
PRF = Ba;                   % 脉冲重复频率（Hz），设置为等于多普勒带宽以满足奈奎斯特采样定理
PRT = 1/PRF;                % 脉冲重复时间（s）
ds = PRT;                   % 慢时间域采样间隔（s），等于脉冲重复时间
Nslow = ceil((Xmax-Xmin+Lsar)/V/ds);  % 慢时间域所需采样点数（覆盖目标区域+合成孔径长度）
Nslow = 2^nextpow2(Nslow);             % 将采样点数调整为2的幂次方，以便FFT高效计算
sn = linspace((Xmin-Lsar/2)/V, (Xmax+Lsar/2)/V, Nslow);  % 慢时间离散序列（s），覆盖整个合成孔径
PRT = (Xmax-Xmin+Lsar)/V/Nslow;       % 刷新脉冲重复时间（s）
PRF = 1/PRT;                           % 刷新脉冲重复频率（Hz）
ds = PRT;                              % 刷新慢时间采样间隔（s）

%%================================================================
%% 快时间域（距离向）参数设计
Tr = 5e-6;                  % 发射脉冲持续时间，5 μs
Br = 30e6;                  % 线性调频（LFM）信号带宽，30 MHz
Kr = Br/Tr;                 % 调频斜率（Hz/s），即LFM的频率变化速率
Fsr = 3*Br;                 % 距离向采样频率（Hz），取带宽的3倍以防止混叠
dt = 1/Fsr;                 % 距离向采样间隔（s）
Rmin = sqrt((Yc-Y0)^2 + H^2);                % 场景近端最小斜距（m）
Rmax = sqrt((Yc+Y0)^2 + H^2 + (Lsar/2)^2);  % 场景远端最大斜距（含方位偏移）（m）
Nfast = ceil(2*(Rmax-Rmin)/C/dt + Tr/dt);    % 距离向采样点数（覆盖双程时延差及脉冲宽度）
Nfast = 2^nextpow2(Nfast);                    % 调整为2的幂次方
tm = linspace(2*Rmin/C, 2*Rmax/C+Tr, Nfast); % 快时间离散序列（s），表示回波双程时延
dt = (2*Rmax/C + Tr - 2*Rmin/C)/Nfast;       % 刷新距离向采样间隔（s）
Fsr = 1/dt;                                   % 刷新距离向采样频率（Hz）

%%================================================================
%% 分辨率计算
DY = C/2/Br;                % 距离分辨率（m），由信号带宽决定
DX = D/2;                   % 方位分辨率（m），等于天线长度的一半

%%================================================================
%% 点目标设置
Ntarget = 2;                % 参与仿真的目标数量
% 目标格式：[方位坐标x(m), 地距y(m), 反射率]
Ptarget = [Xmin,      Yc,          1   % 目标1：方位Xmin，距离Yc（中心）
           Xmin,      Yc+10*DY,    1   % 目标2：同方位，距离Yc+10个分辨单元
           Xmin+20*DX, Yc+50*DY,  1]; % 目标3（备用，当Ntarget≥3时使用）

%% 打印参数信息
disp('Parameters:')
disp('Sampling Rate in fast-time domain');  disp(Fsr/Br)   % 距离向过采样率
disp('Sampling Number in fast-time domain'); disp(Nfast)   % 距离向采样点数
disp('Sampling Rate in slow-time domain');  disp(PRF/Ba)   % 方位向过采样率
disp('Sampling Number in slow-time domain'); disp(Nslow)   % 方位向采样点数
disp('Range Resolution');    disp(DY)        % 距离分辨率（m）
disp('Cross-range Resolution'); disp(DX)    % 方位分辨率（m）
disp('SAR integration length'); disp(Lsar)  % 合成孔径长度（m）
disp('Position of targets');  disp(Ptarget) % 目标坐标

%%================================================================
%% 原始回波信号生成
K = Ntarget;                % 目标总数
N = Nslow;                  % 慢时间采样点数（方位向）
M = Nfast;                  % 快时间采样点数（距离向）
T = Ptarget;                % 目标坐标矩阵
Srnm = zeros(N, M);         % 初始化回波数据矩阵（N行×M列，行=方位，列=距离）

for k = 1:1:K               % 逐个目标叠加回波
    sigma = T(k, 3);                     % 第k个目标的后向散射系数（反射率）
    Dslow = sn*V - T(k, 1);             % 方位向位移：平台位置 - 目标方位坐标（m）
    R = sqrt(Dslow.^2 + T(k,2)^2 + H^2); % 每个慢时间采样时刻到目标的瞬时斜距（m）
    tau = 2*R/C;                         % 每个慢时间对应的回波双程时延（s）
    Dfast = ones(N,1)*tm - tau'*ones(1,M); % 快时间矩阵：各行为(tm - τ_k(sn))（s）
    phase = pi*Kr*Dfast.^2 - (4*pi/lambda)*(R'*ones(1,M)); % 回波相位：LFM相位 - 距离相位
    % 回波信号：在脉冲持续时间内（0<Dfast<Tr）且在合成孔径范围内叠加
    Srnm = Srnm + sigma * exp(j*phase) .* (0<Dfast & Dfast<Tr) ...
           .* ((abs(Dslow)<Lsar/2)' * ones(1,M));
end

%%================================================================
%% 距离压缩（Range Compression）
tr = tm - 2*Rmin/C;         % 距离压缩参考时间轴（消去最小双程时延偏置）
Refr = exp(j*pi*Kr*tr.^2) .* (0<tr & tr<Tr); % 距离向匹配滤波器（参考信号，LFM）
% 频域匹配滤波：回波行FFT × 参考信号共轭FFT，再IFFT得到距离压缩结果
Sr = ifty(fty(Srnm) .* (ones(N,1) * conj(fty(Refr))));
Gr = abs(Sr);               % 距离压缩后信号的幅度

%% 方位压缩（Azimuth Compression）
ta = sn - Xmin/V;           % 方位向参考时间轴（以成像区域起始位置为基准）
Refa = exp(j*pi*Ka*ta.^2) .* (abs(ta) < Tsar/2); % 方位向匹配滤波器（参考信号）
% 频域匹配滤波：回波列FFT × 参考信号共轭FFT（转置后扩展为矩阵），再IFFT
Sa = iftx(ftx(Sr) .* (conj(ftx(Refa)).' * ones(1,M)));
Ga = abs(Sa);               % 方位压缩后信号的幅度（最终成像结果）

%%================================================================
%% 绘制成像强度图
colormap(gray);             % 使用灰度色图
figure(1)
subplot(211);
row = tm*C/2 - 2008;        % 距离坐标轴（m）
col = sn*V - 26;            % 方位坐标轴（m）
imagesc(row, col, 255-Gr);  % 距离压缩后图像（翻转灰度使目标显示为亮点）
axis([Yc-Y0, Yc+Y0, Xmin-Lsar/2, Xmax+Lsar/2]);   % 设置坐标轴范围
xlabel('\rightarrow\itRange in meters');            % 横轴：距离（m）
ylabel('\itAzimuth in meters\leftarrow');           % 纵轴：方位（m）
title('Stripmap SAR after range compression');      % 标题：距离压缩后

subplot(212);
imagesc(row, col, 255-Ga);  % 方位压缩后图像（最终成像结果）
axis([Yc-Y0, Yc+Y0, Xmin-Lsar/2, Xmax+Lsar/2]);
xlabel('\rightarrow\itRange in meters');
ylabel('\itAzimuth in meters\leftarrow');
title('Stripmap SAR after range and azimuth compression'); % 标题：成像结果

%%================================================================
%% 绘制三维瀑布图
figure(2)
waterfall(real(Srnm(200:205, :)));  % 原始回波信号（第200~205行）的实部瀑布图
axis tight;
xlabel('Range');            % 距离轴
ylabel('Azimuth');          % 方位轴
title('Real part of the raw signal'); % 原始回波信号实部

figure(3)
waterfall(Gr(200:205, 600:1000));   % 距离压缩后信号（局部区域）的瀑布图
axis tight;
xlabel('Range');
ylabel('Azimuth');
title('Stripmap SAR after range compression'); % 距离压缩后

figure(4)
mesh(Ga(200:300, 750:860));         % 方位压缩后信号（局部区域）的网格图
axis tight;
xlabel('Range');
ylabel('Azimuth');
title('Stripmap SAR after range and azimuth compression'); % 成像结果

%%================================================================
%% 绘制-3dB分辨率轮廓图
figure(5)
a = max(max(Ga));           % 成像结果的最大幅度值
% 绘制幅度为最大值的0.707（即-3dB）和1.0（最大值）两条等高线
contour(row, col, Ga, [0.707*a, a], 'b');
grid on;
axis([9995, 10050, -20, 20]);       % 限制显示范围，聚焦目标区域
xlabel('\rightarrow\itRange in meters');
ylabel('\itAzimuth in meters\leftarrow');
title('Resolution Demo: -3dB contour');  % 标题：-3dB分辨率演示
%%================================================================
