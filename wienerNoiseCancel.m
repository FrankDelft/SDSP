clear;close all;

[clean,fs] = audioread("clean_speech.wav");
[noise,fs2] = audioread("Speech_shaped_noise.wav");
[babble_noise,fs3] = audioread("babble_noise.wav");
white_noise = wgn(633000,1,0);
t = [1:633000]';
non_stationary_noise = 1/2*t/200000 + (t/200000).*randn(size(t));
compound_noise = noise(1:633000)+white_noise*0.05+0.2*sin(t)+0.3*sin(3.*t);

%construct a clean signal of the same length for result generation
clean_speech=compose_signal(clean,zeros(633000,1));

% construct the noisy signals
noise_speech_signal=compose_signal(clean,noise);
babble_speech_signal=compose_signal(clean,babble_noise);
white_noise_speech_signal=compose_signal(clean,white_noise);
non_stationary_noise_speech_signal=compose_signal(clean,non_stationary_noise);

L=600;

% plot clean speech
subplot(6,3,[1,2,3])
plot(clean_speech)

% plot noises
subplot(6,3,4)
plot(noise)
subplot(6,3,7)
plot(babble_noise)
subplot(6,3,10)
plot(white_noise)
subplot(6,3,13)
plot(non_stationary_noise)
subplot(6,3,16)
plot(compound_noise)

% plot noisy signal
subplot(6,3,5)
plot(noise_speech_signal)
subplot(6,3,8)
plot(babble_speech_signal)
subplot(6,3,11)
plot(white_noise_speech_signal)
subplot(6,3,14)
plot(non_stationary_noise_speech_signal)
subplot(6,3,17)
plot(noise_speech_signal)

% Perform noise cancelation and print metrics
disp('Speech shaped noise results before filtering')
metrics(noise_speech_signal, clean_speech);
n_c=noiseCancel(noise_speech_signal, noise, L, 0.0);
disp('Speech shaped noise results after filtering')
metrics(n_c, clean_speech);
subplot(6,3,6)
plot(n_c)

disp('Babble noise results before filtering')
metrics(babble_speech_signal, clean_speech);
n_c=noiseCancel(babble_speech_signal, babble_noise, L, 0.0);
disp('Babble noise results after filtering')
metrics(n_c, clean_speech);
subplot(6,3,9)
plot(n_c)

disp('white noise results before filtering')
metrics(white_noise_speech_signal, clean_speech);
n_c=noiseCancel(white_noise_speech_signal, white_noise, L, 0.0);
disp('white noise results after filtering')
metrics(n_c, clean_speech);
subplot(6,3,12)
plot(n_c)

disp('non stationary noise results before filtering')
metrics(non_stationary_noise_speech_signal, clean_speech);
n_c=noiseCancel(non_stationary_noise_speech_signal, non_stationary_noise, L, 0.0);
disp('non stationary noise results after filtering')
metrics(n_c, clean_speech);
subplot(6,3,15)
plot(n_c)

disp('compound noise results before filtering')
metrics(noise_speech_signal, clean_speech);
n_c=noiseCancel(noise_speech_signal, compound_noise, L, 0.0);
disp('compound noise results after filtering')
metrics(n_c, clean_speech);
subplot(6,3,18)
plot(n_c)

% soundsc(n_c,fs)

function noiseCancelled = noiseCancel(audioData, noiseData, L, Overlap)
    % Calculate number of segments
    N=length(audioData);
    D=floor(L*(1-Overlap));
    K=floor((N-L+D)/D);
    % N_new=K*D-D+L;
    
    % Initialize segmented data
    segmented_data = zeros(L, K);
    segmented_noise = zeros(L, K);
    
    filtered_signal = zeros(size(audioData));
    % Segment audio data and compute periodograms
    for i = 1:K
        % Extract segment
        frame_start=1 + (i - 1) * (D);
        frame_end=(i-1) * (D)+L;
        segmented_data(:, i) = audioData(frame_start : frame_end);
        segmented_noise(:, i) = noiseData(frame_start : frame_end);
        r_n = xcorr(segmented_noise(:, i), 'biased');
        r_ny = xcorr(segmented_noise(:, i), segmented_data(:, i), 'biased');

        R_n_matrix = toeplitz(r_n(L:end));
        w = R_n_matrix \ r_ny(L:end);

        estimated_noise = conv(segmented_noise(:, i), w, 'full');
        audio_canceled=segmented_data(:, i)-estimated_noise(1:L);
        filtered_signal(frame_start:frame_end) = filtered_signal(frame_start:frame_end) + audio_canceled;
    end
    noiseCancelled= filtered_signal;
end

function noisey_signal = compose_signal(clean, noise)
    noisey_signal=clean+noise(1:length(clean));
end