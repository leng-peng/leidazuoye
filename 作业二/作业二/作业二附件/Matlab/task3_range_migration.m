%%================================================================
%% 文件名: task3_range_migration.m
%% 作业三：单点目标 —— 距离弯曲（Range Cell Migration, RCM）分析
%%
%% 内容：
%%   1. 计算理论距离弯曲量及其影响因素
%%   2. 在距离压缩后的数据中可视化距离弯曲现象
%%   3. 实现距离弯曲校正（RCMC），对比校正前后的成像结果
%%   4. 讨论速度、距离、带宽等参数对距离弯曲的影响
%%================================================================
clear;clc;close all;

%%================================================================
%% 系统参数
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
%% 理论距离弯曲量分析
%%
%% 距离弯曲原理：
%%   SAR平台飞行时，目标到平台的斜距 R(s) 随慢时间 s 变化：
%%     R(s) = sqrt( (V*s - x0)^2 + y0^2 + H^2 )
%%   在最近点（s = x0/V）时斜距最小，记为 R_closest
%%   距离弯曲量定义为：
%%     ΔR(s) = R(s) - R_closest
%%   合成孔径两端（s = ±Tsar/2）处最大距离弯曲：
%%     ΔR_max = R(Tsar/2) - R_closest ≈ (V*Tsar/2)^2 / (2*R_closest)
%%            = (Lsar/2)^2 / (2*R0)
%%            = lambda^2 * R0 / (8 * D^2)
%%
%% 当 ΔR_max > DY（距离分辨率）时，距离弯曲将跨越多个距离单元，
%% 必须进行距离弯曲校正（RCMC）才能得到聚焦良好的图像。

DRmax_theory = (Lsar/2)^2 / (2*R0);  % 理论最大距离弯曲量（m）
DRcells = DRmax_theory / DY;          % 弯曲量对应的距离分辨单元数

disp('=== 作业三：距离弯曲分析 ===');
disp(['合成孔径长度 Lsar       = ', num2str(Lsar,'%.2f'), ' m']);
disp(['最短斜距 R0             = ', num2str(R0,'%.2f'), ' m']);
disp(['距离分辨率 DY           = ', num2str(DY,'%.2f'), ' m']);
disp(['理论最大距离弯曲量      = ', num2str(DRmax_theory,'%.2f'), ' m']);
disp(['对应距离分辨单元数      = ', num2str(DRcells,'%.2f'), ' cells']);
if DRcells > 1
    disp('>> 距离弯曲 > 1个分辨单元，需进行RCMC！');
else
    disp('>> 距离弯曲 < 1个分辨单元，可忽略。');
end

%% 影响因素分析：改变各参数观察 ΔR_max 的变化
disp(' ');
disp('--- 距离弯曲影响因素分析 ---');
for Vi = [50, 100, 200]
    R0i   = R0;
    Lsari = lambda*R0i/D;
    DRi   = (Lsari/2)^2 / (2*R0i);
    fprintf('V = %4d m/s, Lsar = %.1f m, ΔRmax = %.2f m (%.2f cells)\n', ...
            Vi, lambda*R0i/D*(Vi/V), DRi*(Vi/V)^2, DRi*(Vi/V)^2/DY);
end
for Yci = [5000, 10000, 20000]
    R0i   = sqrt(Yci^2 + H^2);
    Lsari = lambda*R0i/D;
    DRi   = (Lsari/2)^2 / (2*R0i);
    fprintf('Yc = %5d m, R0 = %.1f m, ΔRmax = %.2f m (%.2f cells)\n', ...
            Yci, R0i, DRi, DRi/DY);
end

%%================================================================
%% 单点目标仿真
Ntarget = 1;
Ptarget = [Xmin, Yc, 1];   % 单个点目标位于成像区域中心

%% 回波生成
N    = Nslow;
M    = Nfast;
Srnm = zeros(N, M);
sigma = Ptarget(1,3);
Dslow = sn*V - Ptarget(1,1);
R_inst = sqrt(Dslow.^2 + Ptarget(1,2)^2 + H^2);  % 各慢时刻的瞬时斜距
tau   = 2*R_inst/C;
Dfast = ones(N,1)*tm - tau'*ones(1,M);
phase = pi*Kr*Dfast.^2 - (4*pi/lambda)*(R_inst'*ones(1,M));
Srnm  = sigma*exp(j*phase) .* (0<Dfast & Dfast<Tr) ...
        .* ((abs(Dslow)<Lsar/2)'*ones(1,M));

%%================================================================
%% 距离压缩（无RCMC）
tr   = tm - 2*Rmin/C;
Refr = exp(j*pi*Kr*tr.^2) .* (0<tr & tr<Tr);
Sr   = ifty(fty(Srnm) .* (ones(N,1)*conj(fty(Refr))));
Gr   = abs(Sr);  % 距离压缩后幅度

%%================================================================
%% RCMC（距离弯曲校正）
%%
%% RD域RCMC原理：
%%   1. 对距离压缩后的数据做方位向FFT → 进入距离-多普勒域 (ft, fa)
%%   2. 对每个多普勒频率 fa，距离弯曲量为：
%%        ΔR(fa) = R0 / sqrt(1 - (lambda*fa/(2V))^2) - R0
%%      对应双程时延差：Δτ(fa) = 2*ΔR(fa)/C
%%      对应距离向采样点偏移：ΔN(fa) = Δτ(fa)/dt
%%   3. 使用线性插值在距离向对每列进行移位校正
%%   4. 方位向IFFT返回距离-慢时间域，再进行方位压缩

%% 方位向FFT → 距离-多普勒域
Sr_rd = fftshift(fft(Sr, [], 1), 1);  % 方位向FFT，以0频为中心

%% 多普勒频率轴
fa_axis = linspace(-PRF/2, PRF/2, N);  % 多普勒频率轴 (Hz)

%% 逐多普勒频率进行距离弯曲校正（线性插值）
Sr_rcmc = zeros(N, M);
fast_time_axis = 0:M-1;    % 快时间索引轴

for i = 1:N
    fa  = fa_axis(i);
    tmp = 1 - (lambda*fa/(2*V))^2;
    if tmp <= 0
        Sr_rcmc(i,:) = Sr_rd(i,:);  % 超出多普勒带宽，不处理
        continue;
    end
    DeltaR = R0/sqrt(tmp) - R0;          % 当前多普勒频率对应的距离弯曲量 (m)
    DeltaN = 2*DeltaR/C/dt;              % 对应的快时间采样点偏移量
    % 线性插值：将每行向距离更小方向移位 DeltaN 个采样点
    new_axis = fast_time_axis + DeltaN;  % 修正后的快时间索引
    % 对有效范围内的数据进行插值
    Sr_rcmc(i,:) = interp1(fast_time_axis, Sr_rd(i,:), new_axis, 'linear', 0);
end

%% 方位向IFFT → 回到距离-慢时间域
Sr_corrected = ifft(ifftshift(Sr_rcmc, 1), [], 1);

%% 方位压缩（无RCMC）
ta   = sn - Xmin/V;
Refa = exp(j*pi*Ka*ta.^2) .* (abs(ta)<Tsar/2);
Sa_noRCMC = iftx(ftx(Sr) .* (conj(ftx(Refa)).'*ones(1,M)));
Ga_noRCMC = abs(Sa_noRCMC);

%% 方位压缩（有RCMC）
Sa_RCMC = iftx(ftx(Sr_corrected) .* (conj(ftx(Refa)).'*ones(1,M)));
Ga_RCMC = abs(Sa_RCMC);

%%================================================================
%% 坐标轴定义
row = tm*C/2;   % 斜距坐标（m）
col = sn*V;     % 方位坐标（m）

%%================================================================
%% 图1：原始回波信号（部分行），观察距离弯曲曲线
figure(1)
% 选择目标方位位置附近的行
az_center = round(Nslow/2);
az_range  = az_center + (-60:60);
imagesc(row, col(az_range), Gr(az_range,:));
colormap(hot); colorbar;
xlabel('斜距 (m)'); ylabel('方位 (m)');
title('距离压缩后图像（显示距离弯曲曲线，无RCMC）');
% 叠加理论距离弯曲曲线
R_closest = sqrt(Ptarget(1,2)^2 + H^2);
s_arr = (col(az_range) - Ptarget(1,1)) / V;
R_curve = sqrt((s_arr*V - Ptarget(1,1)).^2 + Ptarget(1,2)^2 + H^2);
hold on;
plot(R_curve, col(az_range), 'g--', 'LineWidth', 2, 'DisplayName', '理论距离弯曲曲线');
legend('Location','best');

%%================================================================
%% 图2：距离-多普勒域中的距离弯曲（RCMC前后对比）
figure(2)
fa_axis_plot = linspace(-PRF/2, PRF/2, N);

subplot(121)
imagesc(row, fa_axis_plot, abs(Sr_rd));
colormap(hot); colorbar;
xlabel('斜距 (m)'); ylabel('多普勒频率 (Hz)');
title('距离-多普勒域（RCMC前）');
axis([R_closest-100, R_closest+100, -PRF/4, PRF/4]);

subplot(122)
imagesc(row, fa_axis_plot, abs(Sr_rcmc));
colormap(hot); colorbar;
xlabel('斜距 (m)'); ylabel('多普勒频率 (Hz)');
title('距离-多普勒域（RCMC后）');
axis([R_closest-100, R_closest+100, -PRF/4, PRF/4]);

%%================================================================
%% 图3：成像结果对比（无RCMC vs 有RCMC）
colormap(gray);
figure(3)
subplot(121)
imagesc(row, col, 255-Ga_noRCMC);
axis([R_closest-50, R_closest+50, -Lsar/2, Lsar/2]);
xlabel('斜距 (m)'); ylabel('方位 (m)');
title('成像结果（无RCMC）');

subplot(122)
imagesc(row, col, 255-Ga_RCMC);
axis([R_closest-50, R_closest+50, -Lsar/2, Lsar/2]);
xlabel('斜距 (m)'); ylabel('方位 (m)');
title('成像结果（有RCMC）');

%%================================================================
%% 图4：距离向和方位向剖面对比
figure(4)
subplot(211)
% 方位向剖面（在目标距离处截取）
[~, ridx] = min(abs(row - R_closest));
az_profile_noRCMC = Ga_noRCMC(:, ridx);
az_profile_RCMC   = Ga_RCMC(:, ridx);
plot(col, az_profile_noRCMC/max(az_profile_noRCMC), 'b-', 'LineWidth', 1.5, 'DisplayName', '无RCMC');
hold on;
plot(col, az_profile_RCMC/max(az_profile_RCMC), 'r--', 'LineWidth', 1.5, 'DisplayName', '有RCMC');
axis([-Lsar/2, Lsar/2, 0, 1.1]);
xlabel('方位 (m)'); ylabel('归一化幅度');
title('方位向剖面对比（无RCMC vs 有RCMC）');
legend; grid on;

subplot(212)
% 距离向剖面（在目标方位处截取）
[~, aidx] = min(abs(col - Ptarget(1,1)));
rng_profile_noRCMC = Ga_noRCMC(aidx, :);
rng_profile_RCMC   = Ga_RCMC(aidx, :);
plot(row, rng_profile_noRCMC/max(rng_profile_noRCMC), 'b-', 'LineWidth', 1.5, 'DisplayName', '无RCMC');
hold on;
plot(row, rng_profile_RCMC/max(rng_profile_RCMC), 'r--', 'LineWidth', 1.5, 'DisplayName', '有RCMC');
axis([R_closest-50, R_closest+50, 0, 1.1]);
xlabel('斜距 (m)'); ylabel('归一化幅度');
title('距离向剖面对比（无RCMC vs 有RCMC）');
legend; grid on;

%%================================================================
%% 图5：-3dB轮廓线对比
figure(5)
a1 = max(max(Ga_noRCMC));
a2 = max(max(Ga_RCMC));
subplot(121)
contour(row, col, Ga_noRCMC, [0.707*a1, a1], 'b'); grid on;
axis([R_closest-20, R_closest+20, -20, 20]);
xlabel('斜距 (m)'); ylabel('方位 (m)');
title('-3dB等高线（无RCMC）');

subplot(122)
contour(row, col, Ga_RCMC, [0.707*a2, a2], 'r'); grid on;
axis([R_closest-20, R_closest+20, -20, 20]);
xlabel('斜距 (m)'); ylabel('方位 (m)');
title('-3dB等高线（有RCMC）');

%%================================================================
%% 打印距离弯曲影响分析结论
disp(' ');
disp('=== 距离弯曲对成像的影响分析 ===');
disp(['理论距离弯曲量 ΔRmax = ', num2str(DRmax_theory,'%.2f'), ...
      ' m = ', num2str(DRcells,'%.2f'), ' 个距离分辨单元']);
disp('影响：');
disp('  1. 不做RCMC时，同一目标的回波跨越多个距离单元，');
disp('     导致方位压缩匹配滤波失配，方位向分辨率下降、旁瓣升高。');
disp('  2. RCMC将各多普勒频率的信号对齐到同一距离单元，');
disp('     恢复匹配滤波的相干性，显著改善方位聚焦质量。');
disp('影响距离弯曲量的主要因素：');
disp('  - 平台速度 V：ΔRmax ∝ V^2（速度越大，弯曲越严重）');
disp('  - 目标距离 R0：ΔRmax ∝ R0（距离越远，弯曲越严重）');
disp('  - 天线长度 D：ΔRmax ∝ 1/D^2（天线越短，孔径越大，弯曲越严重）');
disp('  - 载波频率 Fc：ΔRmax ∝ lambda^2 ∝ 1/Fc^2（频率越高，弯曲越轻）');
%%================================================================
