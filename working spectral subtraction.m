clear;close all;clc;
[clean,fs] = audioread("clean_speech.wav");
[noise,fs2] = audioread("Speech_shaped_noise.wav");
%[noise,fs2] = audioread("babble_noise.wav");

%variables defined to construct noisy speech with period of noise before
%speech starts
length_noise_seconds=3.5;
N_noise=fs*length_noise_seconds;
N_clean=500000;
%construct noisy speech
noisy_speech=zeros(N_noise+N_clean,1);
noisy_speech(1:N_noise)=noise(1:N_noise);
noisy_speech(N_noise+1:N_clean+N_noise)=clean(1:N_clean)+noise(N_noise+1:N_clean+N_noise);

L=600;

win=ones(L,1);

SNR=snr(noisy_speech, noise(1:N_clean+N_noise));
[spectral_subtract_enhanced,psd_n_est]=spectral_subtraction(noisy_speech,N_noise,L,fs,0.5,win);
SNR2=snr(spectral_subtract_enhanced, noise(1:N_clean+N_noise)); 
audiowrite("sdspProj\signals\spectral sub\spec1.wav",real(spectral_subtract_enhanced),fs);



% Create a figure
figure;

% Create the first subplot
subplot(2, 2, 1);  % 2 rows, 2 columns, first subplot
plot(real(spectral_subtract_enhanced));
title('enhanced');

% Create the second subplot
subplot(2, 2, 2);  % 2 rows, 2 columns, second subplot
plot(noisy_speech);
title('speech with noise');

% Create the third subplot
subplot(2, 2, 3);  % 2 rows, 2 columns, third subplot
plot(clean(1:N_clean));
title('Clean speech');

% Create the fourth subplot
ifft(sqrt(psd_n_est*L))
subplot(2, 2, 4);  % 2 rows, 2 columns, fourth subplot
plot(psd_n_est);
title('Noise removed');
%[freq,PSD_N]=welsch(noise, L, fs,0.5,win);
%plot(linspace(0,1,L/2),10*log10(PSD_N(1:L/2)));

%hold on
%plot(linspace(0,1,L/2),10*log10(psd_n_est(1:L/2)));





function [reconstructed_signal, psd_noise_est] = spectral_subtraction(audioData, noise_length_start, L, fs,Overlap,window)
    % Calculate number of segments
    N=length(audioData);
    D=floor(L*(1-Overlap));
    K=floor((N-L+D)/D);
    N_new=K*D-D+L;
    %calculate energy of the window
    U=sum(window.^2).*1/length(window);
    
    %estimate noise based on starting noise samples
    [~, psd_noise_est]=welsch(audioData(1:noise_length_start),L,fs,0.5,window);
    

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

        alpha=4;
        beta=0.018;
    
        spectral_sub = psd_speech_frame-alpha.*psd_noise_est';
        factor_reduction = beta.*(psd_noise_est');
        %remove musical noise
        enhanced_psd = max(spectral_sub,factor_reduction);

        gain=psd_speech_frame./((psd_noise_est').*100+psd_speech_frame); 
        %enhanced_psd=gain.*psd_speech_frame;

        %now to reconstruct the speech signal
        enhanced_frame=ifft(sqrt(enhanced_psd.*(U*L)).*exp(phases_frame*1i));
        enhancedSpeech(start_frame:end_frame)=enhancedSpeech(start_frame:end_frame)+enhanced_frame./window;

  
    end
    reconstructed_signal=enhancedSpeech;
end




% Periodogram averaging Welsch
function [freq,PSD] = welsch(audioData,L,fs,Overlap,window)
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
    
    % Create a frequency vector
    freq = (0:L - 1) * fs / L;
end

