%% LFM脉冲压缩雷达仿真
% 线性调频（LFM）脉冲压缩雷达仿真函数
% 参数说明：
%   T    - chirp信号的持续脉宽（秒）
%   B    - chirp信号的调频带宽（Hz）
%   Rmin - 观测目标距雷达的最近位置（米）
%   Rmax - 观测目标距雷达的最远位置（米）
%   R    - 一维数组，每个目标相对雷达的斜距（米）
%   RCS  - 一维数组，每个目标的雷达散射截面
%
% 调用示例：
%   LFM_radar(10e-6, 30e6, 10000, 15000, [10500,11000,12000,12005,13000,13002], [1,1,1,1,1,1])

function LFM_radar(T, B, Rmin, Rmax, R, RCS)
if nargin == 0
    T    = 10e-6;                               % 脉冲宽度 10us
    B    = 30e6;                                % 调频带宽 30MHz
    Rmin = 10000; Rmax = 15000;                 % 距离门范围（米）
    R    = [10500, 11000, 12000, 12005, 13000, 13002]; % 目标斜距（米）
    RCS  = [1, 1, 1, 1, 1, 1];                 % 目标雷达散射截面
end

%% 系统参数
C    = 3e8;             % 光速（m/s）
K    = B / T;           % 调频斜率
Rwid = Rmax - Rmin;     % 距离门宽度（米）
Twid = 2 * Rwid / C;    % 接收窗口时宽（秒）
Fs   = 5 * B; Ts = 1 / Fs;          % 采样频率（5倍过采样，保证回波重建精度）
Nwid = ceil(Twid / Ts);              % 接收窗口采样点数

%% 产生回波信号
t  = linspace(2*Rmin/C, 2*Rmax/C, Nwid);   % 接收时间窗口
M  = length(R);                              % 目标个数
td = ones(M, 1) * t - 2*R' / C * ones(1, Nwid);
% 各点目标产生的雷达回波叠加
Srt = RCS * (exp(1j*pi*K*td.^2) .* (abs(td) < T/2));

%% 用FFT/IFFT实现脉冲压缩
Nchirp = ceil(T / Ts);                      % chirp信号采样点数
Nfft   = 2^nextpow2(Nwid + Nchirp - 1);    % FFT点数（线性卷积）

Srw = fft(Srt, Nfft);                       % 回波信号FFT
t0  = linspace(-T/2, T/2, Nchirp);
St  = exp(1j*pi*K*t0.^2);                   % chirp信号复包络
Sw  = fft(St, Nfft);                        % chirp信号FFT
Sot = fftshift(ifft(Srw .* conj(Sw)));      % 脉冲压缩输出

%% 结果显示
N0 = Nfft/2 - Nchirp/2;
Z  = abs(Sot(N0 : N0 + Nwid - 1));
Z  = Z / max(Z);
Z  = 20 * log10(Z + 1e-6);

figure;
subplot(2, 1, 1);
plot(t * 1e6, real(Srt));
axis tight;
xlabel('Time (us)');
ylabel('Amplitude');
title('Radar Echo without Pulse Compression (实部)');
grid on;

subplot(2, 1, 2);
plot(t * C / 2, Z);
axis([Rmin, Rmax, -60, 0]);
xlabel('Range (m)');
ylabel('Amplitude (dB)');
title('Radar Echo after Pulse Compression');
grid on;
end
