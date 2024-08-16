% Define parameters
snr_range = -40:1:20; % SNR range in dB
num_trials = 5; % Number of trials试验次数
fs = 1e3; % Sampling frequency
num_samples = 500; % Number of samples per trial
freq_range = [3, 300]; % Frequency range
signal_frequency = 50; % Example signal frequency

% Define the signals
t = (0:num_samples-1)/fs;
signal = cos(2*pi*signal_frequency*t); % Example signal
noise = randn(1,length(signal));
noise_power = mean(noise.^2);
%figure;
%plot(t, signal, '-o');
%[caf, ~] = cycloacf(noise);
%caf_noise_power = sum(abs(caf).^2);
cov_matrix = calculate_covariance_matrix(noise);
eigenvalues = eig(cov_matrix);
ratio = max(eigenvalues)/min(eigenvalues);

% Placeholder for results
results = struct();

% Loop through each SNR value
for snr_idx = 1:length(snr_range)
    snr = snr_range(snr_idx);
    results(snr_idx).snr = snr;
    
    % Loop through each trial
    
    for trial = 1:num_trials
        % Generate noisy signal
        noisy_signal = add_awgn(signal,snr,noise_power);
        
        % Apply different detection techniques
        results(snr_idx).energy_detection(trial) = energy_detection(noisy_signal);
        results(snr_idx).matched_filter(trial) = matched_filter(noisy_signal, signal);
        results(snr_idx).cyclostationary_detection(trial) = cyclostationary_detection(noisy_signal);
        results(snr_idx).eigenvalue_detection(trial) = eigenvalue_detection(noisy_signal);
    end
    if mod(snr_idx,10) == 0    
        figure;
        plot(t,noisy_signal);
    end
end

% Calculate detection probabilities and false alarm rates
% This is a simplified placeholder, please replace with actual calculation
for snr_idx = 1:length(snr_range)
    results(snr_idx).energy_detection_mean = mean([results(snr_idx).energy_detection]);
    results(snr_idx).matched_filter_mean = mean([results(snr_idx).matched_filter]);
    results(snr_idx).cyclostationary_detection_mean = mean([results(snr_idx).cyclostationary_detection]);
    results(snr_idx).eigenvalue_detection_mean = mean([results(snr_idx).eigenvalue_detection]);
end

% Visualization
figure;
%hold on;
subplot 221;
plot(snr_range, [results.energy_detection_mean], '-o', 'DisplayName', 'Energy Detection');
xlabel('SNR (dB)');
ylabel('Detection Probability');
legend;
title('Detection Probability vs SNR');
grid on;
subplot 222;
plot(snr_range, [results.matched_filter_mean], '-x', 'DisplayName', 'Matched Filter');
xlabel('SNR (dB)');
ylabel('Detection Probability');
legend;
title('Detection Probability vs SNR');
grid on;
subplot 223;
plot(snr_range, [results.cyclostationary_detection_mean], '-s', 'DisplayName', 'Cyclostationary Detection');
xlabel('SNR (dB)');
ylabel('Detection Probability');
legend;
title('Detection Probability vs SNR');
grid on;
subplot 224;
plot(snr_range, [results.eigenvalue_detection_mean], '-d', 'DisplayName', 'Eigenvalue Detection');
xlabel('SNR (dB)');
xlim([-40,20]);
ylabel('Detection Probability');
ylim([0,1]);
legend;
title('Detection Probability vs SNR');
grid on;

%energy_noise = sum(abs(noise).^2);
corre_noise = max(abs(conv(noise,signal)));

% Functions for detection techniques
function pd = energy_detection(signal)
    threshold = 563.7; % Example threshold
    energy = sum(abs(signal).^2);
    pd = energy > threshold;
end

function pd = matched_filter(signal, template)
    threshold = 30; % Example threshold
    correlation = abs(conv(signal, fliplr(template)));
    pd = max(correlation) > threshold;
end

function pd = eigenvalue_detection(signal)
    threshold = 1; % Example threshold
    received_signal = zeros(2,size(signal,2)/2);
    for b = 1:size(signal,2)/2
        for a = 1:2
            received_signal(a,b) = signal(2*a+b-2);
        end
    end
    cov_matrix = calculate_covariance_matrix(received_signal);
    eigenvalues = eig(cov_matrix);
    ratio = max(eigenvalues)/min(eigenvalues);
    % Placeholder for eigenvalue detection
    % Implement actual eigenvalue-based detection
    pd = ratio > threshold; % Random result for placeholder
end

function y = add_awgn(x, snr_db, noise_power)%这个函数可以在恒定噪声功率和恒定信噪比的条件下，加上信号
    % x: 原始信号
    % snr_db: 指定的信噪比，单位为dB
    % noise_power: 固定的噪声功率
    
    % 计算所需的信号功率
    snr_linear = 10^(snr_db / 10);  % 将dB转换为线性值
    signal_power = noise_power * snr_linear;  % 所需的信号功率
    
    % 调整信号的功率
    x_power = mean(x.^2);  % 计算原始信号的功率
    scaling_factor = sqrt(signal_power / x_power);  % 计算缩放因子
    x_scaled = x * scaling_factor;  % 调整后的信号
    
    % 生成高斯白噪声
    noise = sqrt(noise_power) * randn(size(x));
    
    % 生成接收信号
    y = x_scaled + noise;
end

function detected = cyclostationary_detection(signal)
    % CYCLOSTATIONARY_DETECTION 使用循环平稳检测方法判断信号中是否存在信号
    % 输入: 
    %   signal - 输入的信号（包含高斯白噪声）
    % 输出:
    %   detected - 如果检测到信号，返回1；否则返回0

    % 设置检测的阈值，根据实际情况调整
    threshold = 5.540193496142288e+05;

    % 计算循环自相关函数（CAF）
    [caf, ~] = cycloacf(signal);

    % 计算CAF的平方和
    caf_power = sum(abs(caf).^2);

    % 比较CAF的平方和与阈值
    detected = caf_power > threshold;
end

function [caf, lags] = cycloacf(signal)
    % 计算循环自相关函数（CAF）
    
    N = length(signal);
    lags = -N+1:N-1;
    caf = zeros(size(lags));
    
    for k = 1:length(lags)
        lag = lags(k);
        if lag < 0
            caf(k) = sum(signal(1:end+lag).*conj(signal(1-lag:end)));
        else
            caf(k) = sum(signal(1+lag:end).*conj(signal(1:end-lag)));
        end
    end
end

function cov_matrix = calculate_covariance_matrix(received_signal)
    % Calculate the covariance matrix of the received signal
    % Input:
    %   received_signal - a matrix where each row is a signal sample
    % Output:
    %   cov_matrix - the covariance matrix of the received signal
    
    % Ensure the input signal is in the correct form (each row is a signal sample)
    if size(received_signal, 1) < size(received_signal, 2)
        received_signal = received_signal.';
    end
    
    % Subtract the mean from the signal (zero-mean)
    mean_signal = mean(received_signal, 1);
    zero_mean_signal = received_signal - mean_signal;
    
    % Calculate the covariance matrix
    cov_matrix = (zero_mean_signal' * zero_mean_signal) / (size(received_signal, 1) - 1);
end
