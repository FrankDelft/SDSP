[clean,fs] = audioread("clean_speech.wav");
[noise,fs] = audioread("Speech_shaped_noise.wav");

noisy_speech=zeros(633000,1);
noisy_speech(1:56000)=noise(1:56000);
noisy_speech(56001:633000)=clean(1:577000)+noise(56001:633000);


speech_signal = noisy_speech;  
noise_signal = noise(1:633000); 

n_c=noiseCancel(speech_signal, noise_signal, 600,0.5);
%soundsc(n_c,fs)

function noiseCancelled = noiseCancel(audioData, noiseData, L, Overlap)
    % Calculate number of segments
    N=length(audioData);
    D=floor(L*(1-Overlap));
    K=floor((N-L+D)/D);
    N_new=K*D-D+L;
    
    % Initialize segmented data and periodogram matrices
    segmented_data = zeros(L, K);
    segmented_noise = zeros(L, K);
    
    filtered_signal = zeros(size(audioData));
    % Segment audio data and compute periodograms
    for i = 1:2
        % Extract segment
        frame_start=1 + (i - 1) * (D);
        frame_end=(i-1) * (D)+L;
        segmented_data(:, i) = audioData(frame_start : frame_end);
        segmented_noise(:, i) = noiseData(frame_start : frame_end);
        r_n = xcorr(segmented_noise(:, i), 'biased');
        r_ny = xcorr(segmented_noise(:, i), segmented_data(:, i), 'biased');

        R_n_matrix = toeplitz(r_n(L:end));
        w = R_n_matrix \ r_ny(L:end);
        
        
        estimated_noise = conv(segmented_noise(:, i), w, 'same');
        noiseCancelled=estimated_noise;
        audio_sub_v1est=segmented_data(:, i)-estimated_noise;
        filtered_signal(frame_start:frame_end) = filtered_signal(frame_start:frame_end) + audio_sub_v1est;
    end

end