function [RMSE, SNR] = metrics(predicted_signal, true_signal)
    RMSE = rmse(true_signal, predicted_signal);
    SNR = snr(true_signal, predicted_signal-true_signal);
    disp('RMSE:')
    disp(RMSE)
    disp('SNR[dB]')
    disp(SNR)
end