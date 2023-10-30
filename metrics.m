function [RMSE, SNR,STOI_metric] = metrics(predicted_signal, true_signal)
    RMSE = rmse(true_signal, predicted_signal);
    SNR = snr(true_signal, predicted_signal-true_signal);
    STOI_metric= stoi(true_signal, predicted_signal,16000);
    disp('RMSE:')
    disp(RMSE)
    disp('SNR[dB]')
    disp(SNR)
    disp('STOI:')
    disp(STOI_metric)
end