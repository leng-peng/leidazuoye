%% LFM脉冲压缩匹配滤波演示
% 对chirp信号进行匹配滤波，并分析距离分辨率
% 参考：PDF《LFM脉冲压缩雷达仿真》第三、四节

clear; clc; close all;

%% 系统参数
T  = 10e-6;     % 脉冲宽度 10us
B  = 30e6;      % 调频带宽 30MHz
K  = B / T;     % 调频斜率
C  = 3e8;       % 光速（m/s）

Fs = 10 * B; Ts = 1 / Fs;  % 采样频率（过采样以提高精度）
N  = round(T / Ts);
t  = linspace(-T/2, T/2, N);

%% 产生chirp信号和匹配滤波器
St = exp(1j * pi * K * t.^2);      % chirp信号复包络
Ht = exp(-1j * pi * K * t.^2);     % 匹配滤波器（时域翻转共轭）

%% 匹配滤波（时域卷积）
Sot = conv(St, Ht);
L   = 2 * N - 1;
t1  = linspace(-T, T, L);

%% 归一化并转为dB
Z  = abs(Sot); Z = Z / max(Z);
Z  = 20 * log10(Z + 1e-6);

%% 理论sinc函数对照
Z1 = abs(sinc(B .* t1));
Z1 = 20 * log10(Z1 + 1e-6);

%% 以 t*B 归一化时间轴
t1B = t1 * B;

figure('Name', 'LFM Matched Filter Output');

subplot(2, 1, 1);
plot(t1B, Z, 'b', t1B, Z1, 'r--');
axis([-15, 15, -50, inf]);
grid on;
legend('仿真结果', 'sinc理论值');
xlabel('归一化时间  t \times B');
ylabel('幅度 (dB)');
title('LFM信号匹配滤波输出（脉冲压缩）');

subplot(2, 1, 2);   % 局部放大
N0  = round(3 * Fs / B);
idx = (N - N0) : (N + N0);
t2  = B * t1(idx);
plot(t2, Z(idx), 'b', t2, Z1(idx), 'r--');
axis([-inf, inf, -50, inf]);
set(gca, 'YTick', [-13.4, -4, 0], 'XTick', [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]);
% -13.4 dB: first sidelobe level of sinc; -4 dB: 3dB width marker (t = ±1/(2B))
grid on;
legend('仿真结果', 'sinc理论值');
xlabel('归一化时间  t \times B');
ylabel('幅度 (dB)');
title('局部放大');

%% 距离分辨率分析
sigma_R = C / (2 * B);
fprintf('=== 距离分辨率分析 ===\n');
fprintf('调频带宽 B = %.0f MHz\n', B / 1e6);
fprintf('理论距离分辨率 σ_R = C/(2B) = %.2f m\n', sigma_R);
fprintf('压缩比 D = T*B = %.0f\n', T * B);

%% 多目标分辨率仿真
fprintf('\n=== 两点目标分辨能力仿真 ===\n');
delta_R_list = [2, 5, 10, 20];   % 两目标间距（米）
Rmin = 0; Rmax = 100;
Twid = 2 * (Rmax - Rmin) / C;
Nwid = ceil(Twid / Ts);
t_echo = linspace(2*Rmin/C, 2*Rmax/C, Nwid);
Nchirp = ceil(T / Ts);
Nfft   = 2^nextpow2(Nwid + Nchirp - 1);
Sw     = fft(St, Nfft);

figure('Name', '两点目标分辨率');
for k = 1:length(delta_R_list)
    dR  = delta_R_list(k);
    R2  = [50, 50 + dR];   % 两目标位置（米）
    td2 = ones(2, 1) * t_echo - 2 * R2' / C * ones(1, Nwid);
    Srt2 = sum(exp(1j*pi*K*td2.^2) .* (abs(td2) < T/2), 1);

    Srw2 = fft(Srt2, Nfft);
    Sot2 = fftshift(ifft(Srw2 .* conj(Sw)));
    N0 = Nfft/2 - Nchirp/2;
    Zk = abs(Sot2(N0 : N0 + Nwid - 1));
    Zk = Zk / max(Zk);
    Zk = 20 * log10(Zk + 1e-6);

    subplot(2, 2, k);
    plot(t_echo * C / 2, Zk);
    axis([40, 80, -60, 0]);
    xlabel('距离 (m)');
    ylabel('幅度 (dB)');
    title(sprintf('目标间距 = %d m (σ_R = %.0f m)', dR, sigma_R));
    grid on;
    xline(R2(1), 'r--'); xline(R2(2), 'r--');
    legend('压缩输出', '目标真实位置');
end
sgtitle('两点目标距离分辨率仿真');
