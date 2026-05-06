%%================================================================
%% 文件名: task4_same_range_diff_azimuth.m
%% 作业四：设置同距离不同方位多个点目标，进行回波生成及成像仿真。
%%         根据作业三的分析调整距离弯曲量，分别在距离方位域和
%%         距离-多普勒域观察距离弯曲现象。
%%
%% 方案：
%%   通过改变平台速度 V（增大到原来的2倍），使距离弯曲量显著增加
%%   （ΔRmax ∝ V^2，速度加倍则弯曲量增加4倍），以便更清晰地观察。
%%   同时保持标准参数作为对比。
%%================================================================
clear;clc;close all;

%%================================================================
%% 系统参数（基础版本，同原始程序）
C      = 3e8;
Fc     = 1e9;
lambda = C/Fc;
Xmin   = 0;
Xmax   = 50;
Yc     = 10000;
Y0     = 500;
H      = 5000;
D      = 4;
Tr     = 5e-6;
Br     = 30e6;
Kr     = Br/Tr;
DY     = C/2/Br;
DX     = D/2;

%% 定义两种速度场景进行对比
V_list  = [100, 200];       % 原始速度 vs 加倍速度（对应4倍距离弯曲量）
labels  = {'V=100 m/s（原始）', 'V=200 m/s（弯曲量增大4倍）'};

for scene = 1:2
    V = V_list(scene);

    %% 根据速度重新计算参数
    R0   = sqrt(Yc^2 + H^2);
    Lsar = lambda*R0/D;
    Tsar = Lsar/V;
    Ka   = -2*V^2/lambda/R0;
    Ba   = abs(Ka*Tsar);
    PRF  = Ba;
    PRT  = 1/PRF;
    Nslow = ceil((Xmax-Xmin+Lsar)/V/PRT);
    Nslow = 2^nextpow2(Nslow);
    sn    = linspace((Xmin-Lsar/2)/V, (Xmax+Lsar/2)/V, Nslow);
    PRT   = (Xmax-Xmin+Lsar)/V/Nslow;
    PRF   = 1/PRT;
    Fsr   = 3*Br;
    dt    = 1/Fsr;
    Rmin  = sqrt((Yc-Y0)^2 + H^2);
    Rmax  = sqrt((Yc+Y0)^2 + H^2 + (Lsar/2)^2);
    Nfast = ceil(2*(Rmax-Rmin)/C/dt + Tr/dt);
    Nfast = 2^nextpow2(Nfast);
    tm    = linspace(2*Rmin/C, 2*Rmax/C+Tr, Nfast);
    dt    = (2*Rmax/C + Tr - 2*Rmin/C)/Nfast;
    Fsr   = 1/dt;

    %% 理论距离弯曲量
    DRmax = (Lsar/2)^2 / (2*R0);

    fprintf('\n=== %s ===\n', labels{scene});
    fprintf('Lsar = %.1f m, ΔRmax = %.2f m = %.2f cells\n', Lsar, DRmax, DRmax/DY);

    %%================================================================
    %% 目标设置：同距离（Yc），不同方位
    %%   方位间距设为 30 个方位分辨单元，确保可分辨
    Ntarget = 3;
    Ptarget = [ Xmin,           Yc, 1   % 目标1：左方位
                Xmin+15*DX,    Yc, 1   % 目标2：中心方位
                Xmin+30*DX,    Yc, 1]; % 目标3：右方位

    %%================================================================
    %% 回波生成
    N    = Nslow;
    M    = Nfast;
    T    = Ptarget;
    Srnm = zeros(N, M);

    for k = 1:Ntarget
        sigma  = T(k,3);
        Dslow  = sn*V - T(k,1);
        R_inst = sqrt(Dslow.^2 + T(k,2)^2 + H^2);
        tau    = 2*R_inst/C;
        Dfast  = ones(N,1)*tm - tau'*ones(1,M);
        phase  = pi*Kr*Dfast.^2 - (4*pi/lambda)*(R_inst'*ones(1,M));
        Srnm   = Srnm + sigma*exp(j*phase) .* (0<Dfast & Dfast<Tr) ...
                 .* ((abs(Dslow)<Lsar/2)'*ones(1,M));
    end

    %%================================================================
    %% 距离压缩
    tr   = tm - 2*Rmin/C;
    Refr = exp(j*pi*Kr*tr.^2) .* (0<tr & tr<Tr);
    Sr   = ifty(fty(Srnm) .* (ones(N,1)*conj(fty(Refr))));
    Gr   = abs(Sr);

    %% 方位压缩（无RCMC）
    ta   = sn - Xmin/V;
    Refa = exp(j*pi*Ka*ta.^2) .* (abs(ta)<Tsar/2);
    Sa   = iftx(ftx(Sr) .* (conj(ftx(Refa)).'*ones(1,M)));
    Ga   = abs(Sa);

    %% 距离-多普勒域
    Sr_rd = fftshift(fft(Sr, [], 1), 1);  % 距离压缩后，方位FFT → 距离-多普勒域
    fa_axis = linspace(-PRF/2, PRF/2, N); % 多普勒频率轴

    %%================================================================
    %% 坐标轴
    row = tm*C/2;
    col = sn*V;
    R_center = sqrt(Yc^2 + H^2);   % 目标斜距

    %%================================================================
    %% 图(scene*10+1)：距离方位域中的距离弯曲
    colormap(gray);
    figure(scene*10 + 1);
    % 显示距离压缩后的图像
    imagesc(row, col, 255-Gr);
    axis([R_center-200, R_center+200, Xmin-Lsar/2, Xmax+Lsar/2]);
    xlabel('斜距 (m)'); ylabel('方位 (m)');
    title(['距离方位域：距离压缩后（', labels{scene}, '）']);
    colormap(gray);
    hold on;
    % 叠加各目标的理论距离弯曲曲线
    colors_rng = {'g','r','b'};
    for k = 1:Ntarget
        s_arr  = col / V;
        R_curve = sqrt((s_arr*V - T(k,1)).^2 + T(k,2)^2 + H^2);
        % 仅在合成孔径范围内绘制
        mask = abs(s_arr*V - T(k,1)) < Lsar/2;
        plot(R_curve(mask), col(mask), [colors_rng{k},'--'], ...
             'LineWidth', 1.5, 'DisplayName', ['目标',num2str(k),'理论曲线']);
    end
    legend('Location','best');

    %%================================================================
    %% 图(scene*10+2)：距离-多普勒域中的距离弯曲
    figure(scene*10 + 2);
    imagesc(row, fa_axis, abs(Sr_rd));
    axis([R_center-200, R_center+200, -PRF/3, PRF/3]);
    colormap(hot); colorbar;
    xlabel('斜距 (m)'); ylabel('多普勒频率 (Hz)');
    title(['距离-多普勒域：RCMC前（', labels{scene}, '）']);
    % 叠加理论距离弯曲曲线（在RD域中为弧线）
    hold on;
    fa_dense = linspace(-PRF/3, PRF/3, 300);
    for k = 1:Ntarget
        Rk = sqrt(T(k,2)^2 + H^2);
        tmp = max(1 - (lambda*fa_dense/(2*V)).^2, 1e-6);
        R_fa = Rk ./ sqrt(tmp);   % 每个多普勒频率对应的斜距
        plot(R_fa, fa_dense, [colors_rng{k},'--'], 'LineWidth', 1.5, ...
             'DisplayName', ['目标',num2str(k),'理论曲线']);
    end
    legend('Location','best');

    %%================================================================
    %% 图(scene*10+3)：最终成像结果（无RCMC）
    figure(scene*10 + 3);
    colormap(gray);
    imagesc(row, col, 255-Ga);
    axis([R_center-200, R_center+200, Xmin-Lsar/2, Xmax+Lsar/2]);
    xlabel('斜距 (m)'); ylabel('方位 (m)');
    title(['成像结果（无RCMC）——', labels{scene}]);
end

%%================================================================
%% 对比总结：打印两种速度下的距离弯曲量
disp(' ');
disp('=== 同距离不同方位目标——距离弯曲量对比 ===');
for scene = 1:2
    V    = V_list(scene);
    R0   = sqrt(Yc^2 + H^2);
    Lsar = lambda*R0/D;
    DRmax = (Lsar/2)^2 / (2*R0);
    fprintf('%s: ΔRmax = %.2f m = %.2f 个距离分辨单元\n', labels{scene}, DRmax, DRmax/DY);
end
disp(' ');
disp('观察结论：');
disp('  1. 在距离方位域（距离压缩后图像），距离弯曲表现为');
disp('     各目标的能量分布在一条抛物线上，而非水平直线；');
disp('     速度越大、弯曲越明显，成像聚焦越差。');
disp('  2. 在距离-多普勒域，距离弯曲表现为不同多普勒频率');
disp('     对应不同的目标斜距（曲线而非竖直线），速度越大曲率越大。');
disp('  3. 若不进行RCMC，方位压缩时匹配滤波失配，');
disp('     导致方位向旁瓣升高、主瓣展宽、图像散焦。');
%%================================================================
