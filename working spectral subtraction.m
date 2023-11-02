clear;close all;clc;
[clean,fs] = audioread("clean_speech.wav");
[noise,fs2] = audioread("Speech_shaped_noise.wav");
[babble_noise,fs3] = audioread("babble_noise.wav");
clean=clean(1:500000);

white_noise = wgn(633000,1,0);
echo_noise = [zeros(2000,1);0.9.*clean]+0.0001.*wgn(2000+length(clean),1,0);
t = [1:633000]';
non_stationary_noise = 1/2*t/200000 + (t/200000).*randn(size(t));
compound_noise = noise(1:633000)+white_noise*0.05+0.2*sin(t)+0.3*sin(3.*t);

length_noise_seconds=1.2;
N_noise=(fs)*length_noise_seconds;

%construct a clean signal of the same length for result generation
clean_speech=compose_signal(clean,zeros(633000,1),length_noise_seconds,fs);

% construct the noisy signals
noise_speech_signal=compose_signal(clean,noise,length_noise_seconds,fs);
babble_speech_signal=compose_signal(clean,babble_noise,length_noise_seconds,fs);
echo_noise_speech_signal=[clean;zeros(2000,1)]+echo_noise;
non_stationary_noise_speech_signal=compose_signal(clean,non_stationary_noise,length_noise_seconds,fs);

%% choosing window, window length and alpha & beta values

%setup windows and window length
L=640;
%win=ones(L,1);
%win=hann(L);
win=hamming(L);
%win=bartlett(L);



alpha=3;
beta=0.025;


N_clean=500000;
spectral_subtract_enhanced=spectral_subtraction(noise_speech_signal,N_noise,L,0.5,win, alpha,beta);

%audiowrite("bart.wav",spectral_subtract_enhanced,fs);

%metrics(spectral_subtract_enhanced(N_noise+1:end), clean(1:N_clean));
%% Noise estimation
L=640;
win=hamming(L);
psd_true=10*log10(welch(noise,L,0.5,win));

length_noise_seconds=0.6;
N_noise=(fs)*length_noise_seconds;
psd_0_6=10*log10(welch(noise(1:N_noise),L,0.5,win));

length_noise_seconds=1.2;
N_noise=(fs)*length_noise_seconds;
psd_1_2=10*log10(welch(noise(1:N_noise),L,0.5,win));



length_noise_seconds=2.4;
N_noise=(fs)*length_noise_seconds;
psd_2_4=10*log10(welch(noise(1:N_noise),L,0.5,win));

radians=linspace(0,pi,L/2);

psd_est_fig1=figure;
set(psd_est_fig1, 'Position', [50, 50,1800, 400]);
subplot(1,4,1);
plot(radians,psd_0_6(1:L/2))
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel("Frequency(rad/sample)")
ylabel("Magnitude(dB)")
title("Estimate using 0.6s of the noise signal")

subplot(1,4,2);
plot(radians, psd_1_2(1:L/2))
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
title("Estimate using 1.2s of the noise signal")
xlabel("Frequency(rad/sample)")
ylabel("Magnitude(dB)")


subplot(1,4,3);
plot(radians, psd_2_4(1:L/2))
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
title("Estimate using 2.4s of the noise signal")
xlabel("Frequency(rad/sample)")
ylabel("Magnitude(dB)")

subplot(1,4,4);
plot(radians,psd_true(1:L/2))
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
title("Estimate using full noise signal")
xlabel("Frequency(rad/sample)")
ylabel("Magnitude(dB)")

saveas(psd_est_fig1,'./signals/spectral sub/psd_est.png')



%% Now for results
length_noise_seconds=1.2;
N_noise=(fs)*length_noise_seconds;
alpha=3;
beta=0.025;
L=640;
win=hamming(L);

% save compound signals
audiowrite('./signals/spectral sub/speech_shaped_noise_speech_signal.wav', noise_speech_signal, fs)
audiowrite('./signals/spectral sub/babble_speech_signal.wav', babble_speech_signal, fs)
audiowrite('./signals/spectral sub/echo_noise_speech_signal.wav', echo_noise_speech_signal, fs)
audiowrite('./signals/spectral sub/non_stationary_noise_speech_signal.wav', non_stationary_noise_speech_signal, fs)


%setup figure for subplots
fig=figure;
% Set default axes font size
set(0, 'DefaultAxesFontSize', 14);

% Set default title font size
set(0, 'DefaultAxesTitleFontSizeMultiplier', 1.1);
set(fig, 'Position', [50, 50,1200, 1200]);



% plot clean speech
subplot(5,2,[1,2])
plot(clean_speech)
title('Clean Speech');
xlabel('Sample (n)')

% plot noisy signal
subplot(5,2,3)
plot(noise_speech_signal)
title('Clean + Stationary Noise');
xlabel('Sample (n)')
subplot(5,2,5)
plot(babble_speech_signal)
title('Clean + Babble Noise');
xlabel('Sample (n)')
subplot(5,2,7)
plot(echo_noise_speech_signal(50000:end))
title('Clean + Echo Noise');
xlabel('Sample (n)')
subplot(5,2,9)
plot(non_stationary_noise_speech_signal)
title('Clean + Non Stationary Noise');
xlabel('Sample (n)')

scores = zeros(2,10);

% Perform noise cancelation and print metrics
disp('Speech shaped noise results before filtering')
[RMSE, SNR] = metrics(noise_speech_signal, clean_speech);
[scores(1,[1,2])] = [RMSE, SNR];
s_s = spectral_subtraction(noise_speech_signal,N_noise, L, 0.5,win,alpha,beta);
audiowrite('./signals/noise_cancel/noise_cancelation_speech_shaped_noise.wav', s_s, fs);
disp('Speech shaped noise results after filtering')
[RMSE, SNR] = metrics(s_s, clean_speech);
[scores(2,[1,2])] = [RMSE, SNR];
subplot(5,2,4)
plot(s_s)
title('Spectral subtraction of estimated Stationary Noise from noisy signal');
xlabel('Sample (n)')

disp('Babble noise results before filtering')
[RMSE, SNR] = metrics(babble_speech_signal, clean_speech);
[scores(1,[3,4])] = [RMSE, SNR];
s_s = spectral_subtraction(babble_speech_signal,N_noise, L, 0.5,win,alpha,beta);
audiowrite('./signals/noise_cancel/noise_cancelation_babble_noise.wav', s_s, fs);
disp('Babble noise results after filtering')
[RMSE, SNR] = metrics(s_s, clean_speech);
[scores(2,[3,4])] = [RMSE, SNR];
subplot(5,2,6)
plot(s_s)
title('Spectral subtraction of estimated babble noise from noisy signal');
xlabel('Sample (n)')

disp('echo noise results before filtering')
[RMSE, SNR] = metrics(echo_noise_speech_signal(50000:end), [clean(50000:end);zeros(2000,1)]);
[scores(1,[5,6])] = [RMSE, SNR];
s_s = spectral_subtraction(echo_noise_speech_signal(50000:end),N_noise, L, 0.5,win,alpha,beta);
audiowrite('./signals/noise_cancel/noise_cancelation_echo_noise.wav', s_s, fs);
disp('echo noise results after filtering')
[RMSE, SNR] = metrics(s_s, [clean(50000:end);zeros(2000,1)]);
[scores(2,[5,6])] = [RMSE, SNR];
subplot(5,2,8)
plot(s_s)
title('Spectral subtraction of estimated echo noise from noisy signal');
xlabel('Sample (n)')

disp('non-stationary noise results before filtering')
[RMSE, SNR] = metrics(non_stationary_noise_speech_signal, clean_speech);
[scores(1,[7,8])] = [RMSE, SNR];
s_s = spectral_subtraction(non_stationary_noise_speech_signal,N_noise, L, 0.5,win,alpha,beta);
audiowrite('./signals/noise_cancel/noise_cancelation_non_stationary_noise.wav', s_s, fs);
disp('non-stationary noise results after filtering')
[RMSE, SNR] = metrics(s_s, clean_speech);
[scores(2,[7,8])] = [RMSE, SNR];
subplot(5,2,10)
plot(s_s)
title('Spectral subtraction of estimated Non-Stationary noise from noisy signal');
xlabel('Sample (n)')

    subplot(5, 2, [1,2]); % Select the i-th subplot
    xlim([0, 530000]); % Set x-axis limits for the current subplot
for i = 3:10 % Adjust this to the number of subplots in your figure
    subplot(5, 2, i); % Select the i-th subplot
    xlim([0, 530000]); % Set x-axis limits for the current subplot
end

saveas(fig,'./signals/spectral sub/signals.png')


function reconstructed_signal = spectral_subtraction(audioData, noise_length_start, L,Overlap,window, alpha, beta)
    % Calculate number of segments
    N=length(audioData);
    D=floor(L*(1-Overlap));
    K=floor((N-L+D)/D);
    %calculate energy of the window
    U=sum(window.^2).*1/length(window);
    
    %estimate noise based on starting noise samples
    psd_noise_est=welch(audioData(1:noise_length_start),L,0.5,window);
   
    enhancedSpeech=zeros(length(audioData),1);
    % Segment audio data and compute periodograms
    for i = 1:K
        % Extract segment
        start_frame=1 + (i - 1) * (D);
        end_frame=(i-1) * (D)+L;
        frame_data = audioData(start_frame : end_frame);
        % Compute periodogram and store phase information
        fft_segment = fft(frame_data .* window);
        psd_speech_frame = (1/U)*(1 / L) .* abs(fft_segment).^2;
        phases_frame = angle(fft_segment);  % Store phase information

     
    
        spectral_sub = psd_speech_frame-alpha.*psd_noise_est';
        factor_reduction = beta.*(psd_noise_est');
        %reduce musical noise
        enhanced_psd = max(spectral_sub,factor_reduction);

        %now to reconstruct the speech signal
        enhanced_frame=ifft(sqrt(enhanced_psd.*(L)).*exp(phases_frame*1i));
        enhancedSpeech(start_frame:end_frame)=enhancedSpeech(start_frame:end_frame)+enhanced_frame;
    end
    reconstructed_signal=real(enhancedSpeech*Overlap);
end




% Periodogram averaging Welch
function PSD = welch(audioData,L,Overlap,window)
    % percentage overlap
    N=length(audioData);
    D=L*(1-Overlap);
    K=floor((N-L+D)/D);
    segmented_data = zeros(L, K);
    U=sum(window.^2).*1/length(window);

    P_i_hat_w = zeros(L, K);
    for i = 1:K
        segmented_data(:, i) = audioData(1 + (i - 1) * (D) : (i-1) * (D)+L);
        P_i_hat_w(:, i) = (1/U).*(1 / L) .* abs(fft(segmented_data(:, i).*window)).^2;
    end
    
    P_ave_w = zeros(1, L);
    for i = 1:K
        P_ave_w(:) = P_ave_w(:) + P_i_hat_w(:, i);
    end
    PSD = P_ave_w .* (1 / K);

end

function noisy_signal = compose_signal(clean, noise,start_noise_seconds,fs)
    %variables defined to construct noisy speech with period of noise before
    %speech starts
    N_noise=fs*start_noise_seconds;
    N_clean=length(clean);
    
    %construct noisy speech
    noisy_speech=zeros(N_noise+N_clean,1);
    noisy_speech(1:N_noise)=noise(1:N_noise);
    noisy_speech(N_noise+1:N_clean+N_noise)=clean(1:N_clean)+noise(N_noise+1:N_clean+N_noise);
    noisy_signal=noisy_speech;
end

