%%================================================================
%% 文件名: task2_same_azimuth_diff_range.m
%% 作业二：设置同方位不同距离的多个点目标，进行回波生成及成像仿真
%%   - 三个目标：方位坐标均为 Xmin，距离坐标不同
%%   - 给出距离压缩后结果及最终成像结果
%%================================================================
clear;clc;close all;

%%================================================================
%% 系统参数（与 stripmapSAR.m 相同）
C      = 3e8;
Fc     = 1e9;
lambda = C/Fc;
Xmin   = 0;
Xmax   = 50;
Yc     = 10000;
Y0     = 500;
V      = 100;
H      = 5000;
R0     = sqrt(Yc^2 + H^2);
D      = 4;
Lsar   = lambda*R0/D;
Tsar   = Lsar/V;
Ka     = -2*V^2/lambda/R0;
Ba     = abs(Ka*Tsar);
PRF    = Ba;
PRT    = 1/PRF;
ds     = PRT;
Nslow  = ceil((Xmax-Xmin+Lsar)/V/ds);
Nslow  = 2^nextpow2(Nslow);
sn     = linspace((Xmin-Lsar/2)/V, (Xmax+Lsar/2)/V, Nslow);
PRT    = (Xmax-Xmin+Lsar)/V/Nslow;
PRF    = 1/PRT;
ds     = PRT;
Tr     = 5e-6;
Br     = 30e6;
Kr     = Br/Tr;
Fsr    = 3*Br;
dt     = 1/Fsr;
Rmin   = sqrt((Yc-Y0)^2 + H^2);
Rmax   = sqrt((Yc+Y0)^2 + H^2 + (Lsar/2)^2);
Nfast  = ceil(2*(Rmax-Rmin)/C/dt + Tr/dt);
Nfast  = 2^nextpow2(Nfast);
tm     = linspace(2*Rmin/C, 2*Rmax/C+Tr, Nfast);
dt     = (2*Rmax/C + Tr - 2*Rmin/C)/Nfast;
Fsr    = 1/dt;
DY     = C/2/Br;
DX     = D/2;

%%================================================================
%% 目标设置：同方位（x = Xmin = 0），不同距离
%%   目标间距设为 50 个距离分辨单元，确保相互可分辨
Ntarget = 3;
Ptarget = [ Xmin,  Yc - 100*DY,  1   % 目标1：近距端
            Xmin,  Yc,           1   % 目标2：中心
            Xmin,  Yc + 100*DY,  1]; % 目标3：远距端

disp('=== 作业二：同方位不同距离点目标仿真 ===');
disp(['距离分辨率 DY = ', num2str(DY), ' m']);
disp(['方位分辨率 DX = ', num2str(DX), ' m']);
disp('目标位置（方位 x, 地距 y, 反射率）：');
disp(Ptarget);

%%================================================================
%% 回波信号生成
K    = Ntarget;
N    = Nslow;
M    = Nfast;
T    = Ptarget;
Srnm = zeros(N, M);

for k = 1:K
    sigma  = T(k, 3);
    Dslow  = sn*V - T(k, 1);
    R      = sqrt(Dslow.^2 + T(k,2)^2 + H^2);
    tau    = 2*R/C;
    Dfast  = ones(N,1)*tm - tau'*ones(1,M);
    phase  = pi*Kr*Dfast.^2 - (4*pi/lambda)*(R'*ones(1,M));
    Srnm   = Srnm + sigma*exp(j*phase) .* (0<Dfast & Dfast<Tr) ...
             .* ((abs(Dslow)<Lsar/2)'*ones(1,M));
end

%%================================================================
%% 距离压缩
tr   = tm - 2*Rmin/C;
Refr = exp(j*pi*Kr*tr.^2) .* (0<tr & tr<Tr);
Sr   = ifty(fty(Srnm) .* (ones(N,1)*conj(fty(Refr))));
Gr   = abs(Sr);

%% 方位压缩
ta   = sn - Xmin/V;
Refa = exp(j*pi*Ka*ta.^2) .* (abs(ta)<Tsar/2);
Sa   = iftx(ftx(Sr) .* (conj(ftx(Refa)).'*ones(1,M)));
Ga   = abs(Sa);

%%================================================================
%% 坐标轴
row = tm*C/2;   % 斜距坐标（m）
col = sn*V;     % 方位坐标（m）

%% 图1：原始回波信号实部（部分行）
figure(1)
waterfall(real(Srnm(200:205,:)));
axis tight;
xlabel('快时间（距离向）采样点');
ylabel('慢时间（方位向）行号');
title('原始回波信号实部（task2：同方位不同距离三目标）');

%% 图2：距离压缩结果 vs 成像结果
colormap(gray);
figure(2)
subplot(211);
imagesc(row, col, 255-Gr);
axis([Yc-Y0, Yc+Y0, Xmin-Lsar/2, Xmax+Lsar/2]);
xlabel('Range (m)'); ylabel('Azimuth (m)');
title('距离压缩后图像（同方位不同距离三点目标）');

subplot(212);
imagesc(row, col, 255-Ga);
axis([Yc-Y0, Yc+Y0, Xmin-Lsar/2, Xmax+Lsar/2]);
xlabel('Range (m)'); ylabel('Azimuth (m)');
title('成像结果（方位压缩后，无RCMC）');

%% 图3：三目标距离向剖面（验证目标可分辨）
figure(3)
% 取三个目标附近的距离向剖面，选取方位向中间行附近
mid_az = round(Nslow/2);
az_range = mid_az-5 : mid_az+5;
range_profile = max(Ga(az_range, :), [], 1);  % 方位向最大值投影
plot(row, range_profile, 'b-', 'LineWidth', 1.5);
hold on;
% 标注目标理论位置
for k = 1:Ntarget
    Rk = sqrt(T(k,2)^2 + H^2);
    xline(Rk, 'r--', ['目标',num2str(k)], 'LabelOrientation','horizontal');
end
xlabel('斜距 (m)'); ylabel('幅度');
title('距离向剖面（同方位不同距离三目标）');
grid on;

%% 图4：-3dB 分辨率等高线
figure(4)
a = max(max(Ga));
contour(row, col, Ga, [0.707*a, a], 'b');
grid on;
% 显示三个目标区域
Rtargets = arrayfun(@(k) sqrt(T(k,2)^2+H^2), 1:Ntarget);
axis([min(Rtargets)-50, max(Rtargets)+50, -20, 20]);
xlabel('Range (m)'); ylabel('Azimuth (m)');
title('Resolution Demo: -3dB contour（同方位不同距离三目标）');
%%================================================================
