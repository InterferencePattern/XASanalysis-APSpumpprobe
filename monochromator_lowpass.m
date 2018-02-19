function [output_signal] = monochromator_lowpass(input_signal,input_x_values,monochromator_res_eV)
%Monochromator_lowpass Applies a low-pass filter to spectra
%   Blocks high-frequency noise that should be invisible, since it exceeds
%   the monochromator resolution.

transformed_signal = fft(input_signal);
num_bins = length(transformed_signal);

% figure(44);
% plot([0:1/(num_bins/2 -1):1], real(transformed_signal(1:num_bins/2)));

fc = 1/(monochromator_res_eV/1000); % Cutoff frequency in keV
fs = 1/((max(input_x_values)-min(input_x_values))/length(input_x_values)); % Sampling frequency of 1/2 eV

[b,a] = butter(2, fc/(fs/2), 'low');

% H = freqz(b,a, floor(num_bins/2));
% figure(45);
% hold on
% plot([0:1/(num_bins/2 -1):1], abs(H),'r');

output_signal = filter(b,a,input_signal);

figure(46);
plot(input_x_values,output_signal,input_x_values,input_signal);
legend(['Filtered Spectrum','Unfiltered Spectrum'])

end

