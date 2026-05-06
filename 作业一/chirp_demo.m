%% LFM（线性调频）信号仿真演示
% 产生chirp信号，展示其时域波形和幅频特性
% 参考：PDF《LFM脉冲压缩雷达仿真》第二节

clear; clc; close all;

%% 信号参数
T  = 10e-6;     % 脉冲宽度 10us
B  = 30e6;      % 调频带宽 30MHz
K  = B / T;     % 调频斜率

Fs = 2 * B; Ts = 1 / Fs;   % 采样频率和采样间隔
N  = round(T / Ts);         % 采样点数
t  = linspace(-T/2, T/2, N);

%% 产生chirp信号复包络  S(t) = rect(t/T) * exp(j*pi*K*t^2)
St = exp(1j * pi * K * t.^2);

%% 时频分析
freq = linspace(-Fs/2, Fs/2, N);
Sf   = fftshift(fft(St));

%% 绘图
figure('Name', 'LFM Chirp Signal');

subplot(2, 1, 1);
plot(t * 1e6, real(St));
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Real Part of LFM Chirp Signal');
grid on; axis tight;

subplot(2, 1, 2);
plot(freq * 1e-6, abs(Sf));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Magnitude Spectrum of LFM Chirp Signal');
grid on; axis tight;

fprintf('LFM信号参数：\n');
fprintf('  脉冲宽度 T = %.1f us\n', T * 1e6);
fprintf('  调频带宽 B = %.0f MHz\n', B / 1e6);
fprintf('  调频斜率 K = %.3e Hz/s\n', K);
fprintf('  时宽带宽积（压缩比）TB = %.0f\n', T * B);
