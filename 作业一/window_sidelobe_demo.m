%% 脉冲压缩旁瓣抑制：不同窗函数对比
% 研究加权窗函数对匹配滤波旁瓣的抑制效果
% 窗函数：矩形窗、Hamming窗、Hanning窗、Kaiser窗、Chebyshev窗
%
% 参考：作业一第2问 ——
%   研究脉压信号的旁瓣抑制方法，分析不同窗函数对脉压结果的影响

clear; clc; close all;

%% 系统参数
T  = 10e-6;     % 脉冲宽度 10us
B  = 30e6;      % 调频带宽 30MHz
K  = B / T;     % 调频斜率
C  = 3e8;       % 光速

Fs     = 10 * B; Ts = 1 / Fs;
Nchirp = round(T / Ts);         % chirp信号采样点数
t0     = linspace(-T/2, T/2, Nchirp);
St     = exp(1j * pi * K * t0.^2);  % 原始chirp信号

%% 目标场景
Rmin = 10000; Rmax = 15000;
R    = [10500, 11000, 12000, 12005, 13000, 13002];
RCS  = [1, 1, 1, 1, 1, 1];

Rwid = Rmax - Rmin;
Twid = 2 * Rwid / C;
Nwid = ceil(Twid / Ts);
t    = linspace(2*Rmin/C, 2*Rmax/C, Nwid);
M    = length(R);
td   = ones(M, 1) * t - 2 * R' / C * ones(1, Nwid);
Srt  = RCS * (exp(1j*pi*K*td.^2) .* (abs(td) < T/2));  % 回波

Nfft = 2^nextpow2(Nwid + Nchirp - 1);
Srw  = fft(Srt, Nfft);

%% 各窗函数定义
window_list = {
    'Rectangular（矩形窗）',  ones(1, Nchirp);
    'Hamming',               hamming(Nchirp)';
    'Hanning',               hanning(Nchirp)';
    'Kaiser (β=6)',          kaiser(Nchirp, 6)';
    'Chebyshev (60dB)',      chebwin(Nchirp, 60)';
};
nWin = size(window_list, 1);

%% 逐窗函数做脉冲压缩并记录指标
N0 = Nfft/2 - Nchirp/2;
colors = lines(nWin);

figure('Name', '不同窗函数脉冲压缩结果对比', 'Position', [100, 100, 1000, 600]);
ax = axes;
hold on;
peak_sidelobe_dB = zeros(1, nWin);
mainlobe_3dB_m   = zeros(1, nWin);

for k = 1:nWin
    win_name = window_list{k, 1};
    win      = window_list{k, 2};

    % 加窗后的chirp（对参考信号加权）
    St_win = St .* win;

    % 频域匹配滤波
    Sw_win = fft(St_win, Nfft);
    Sot    = fftshift(ifft(Srw .* conj(Sw_win)));

    Z = abs(Sot(N0 : N0 + Nwid - 1));
    Z = Z / max(Z);
    Z_dB = 20 * log10(Z + 1e-6);

    % 峰值旁瓣电平（以最强单点目标Tar1=10500m为参考，排除主瓣±sigma_R区域）
    [~, pk_idx] = max(Z);
    sigma_R = C / (2 * B);
    idx_margin = round(sigma_R / (C * Ts / 2));
    mask = true(1, Nwid);
    mask(max(1, pk_idx - idx_margin) : min(Nwid, pk_idx + idx_margin)) = false;
    if any(mask)
        peak_sidelobe_dB(k) = max(Z_dB(mask));
    end

    % 3dB主瓣宽度（对单点目标）
    half_power = Z_dB(pk_idx) - 3;
    above = Z_dB >= half_power;
    runs  = diff([0, above, 0]);
    start_idx = find(runs == 1);
    end_idx   = find(runs == -1) - 1;
    [~, run_containing] = max(end_idx - start_idx);
    if ~isempty(start_idx)
        width_samp = end_idx(run_containing) - start_idx(run_containing) + 1;
        mainlobe_3dB_m(k) = width_samp * Ts * C / 2;
    end

    plot(ax, t * C / 2, Z_dB, 'Color', colors(k,:), 'LineWidth', 1.2, ...
        'DisplayName', win_name);
end

axis([Rmin, Rmax, -80, 5]);
xlabel('Distance (m)');
ylabel('Amplitude (dB)');
title('脉冲压缩结果：不同窗函数旁瓣抑制效果对比');
legend('Location', 'northeast');
grid on;
hold off;

%% 打印性能指标汇总表
fprintf('\n===== 窗函数对脉压旁瓣的影响 =====\n');
fprintf('%-26s  %s\n', '窗函数', '峰值旁瓣电平 (dB)');
fprintf('%s\n', repmat('-', 1, 50));
for k = 1:nWin
    fprintf('%-26s  %.1f dB\n', window_list{k, 1}, peak_sidelobe_dB(k));
end

%% 单点目标对比：主瓣展宽 vs 旁瓣抑制
figure('Name', '单点目标窗函数对比');
R_single = 12000;
t_s   = linspace(2*Rmin/C, 2*Rmax/C, Nwid);
td_s  = t_s - 2 * R_single / C;
Srt_s = exp(1j*pi*K*td_s.^2) .* (abs(td_s) < T/2);
Srw_s = fft(Srt_s, Nfft);

hold on;
sigma_R = C / (2 * B);
for k = 1:nWin
    win     = window_list{k, 2};
    St_win  = St .* win;
    Sw_win  = fft(St_win, Nfft);
    Sot_s   = fftshift(ifft(Srw_s .* conj(Sw_win)));
    Z_s     = abs(Sot_s(N0 : N0 + Nwid - 1));
    Z_s     = Z_s / max(Z_s);
    Z_s_dB  = 20 * log10(Z_s + 1e-6);
    plot(t_s * C / 2 - R_single, Z_s_dB, 'Color', colors(k, :), ...
        'LineWidth', 1.2, 'DisplayName', window_list{k, 1});
end
axis([-10*sigma_R, 10*sigma_R, -80, 5]);
xlabel('相对距离  Δr (m)');
ylabel('幅度 (dB)');
title(sprintf('单点目标脉压输出（σ_R = %.1f m）', sigma_R));
legend('Location', 'northeast');
grid on;
hold off;

fprintf('\n===== 结论 =====\n');
fprintf('1. 矩形窗：主瓣最窄（最优距离分辨率），旁瓣最高（约-13.4dB）\n');
fprintf('2. Hamming窗：旁瓣约-43dB，主瓣展宽约1.5倍\n');
fprintf('3. Hanning窗：旁瓣约-32dB，主瓣展宽约1.6倍\n');
fprintf('4. Kaiser窗 (β=6)：旁瓣约-44dB，可通过β参数灵活调节旁瓣/分辨率折中\n');
fprintf('5. Chebyshev窗：各旁瓣等幅（-60dB），主瓣展宽最大\n');
fprintf('\n距离分辨率（σ_R = C/(2B) = %.1f m）\n', C/(2*B));
