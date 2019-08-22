clear;
close all;
clc;

% Initialize base parameters
speed_of_light_vacuum = 3e8; % Speed in meters / second
initial_width_constant = 4;
width_step = 0.0001;
supergauss_sigma = 0.3;
supergauss_pow = 8;
threshold = 0.01; % Desired threshold for window search algorithm
SG_order = 2;
SG_window_size = 101;
offset = 0.2;
remove_offset = true;
number_of_steps = 10000;
fwhm_pulse_duration = 20e-15; % Time in seconds
carrier_wavelength = 810e-9; % Wavelength in meters
carrier_frequency = speed_of_light_vacuum / carrier_wavelength;
carrier_angular_frequency = 2 * pi * carrier_frequency;
snr = 20;
number_of_samples = 2048;
cutoff_wavelength = [365 2325] * 1e-9;
cutoff_angular_frequency = wrev(2 * pi * speed_of_light_vacuum ./ cutoff_wavelength);
        
% Generate spectrum
chi = 1/1230;
T = fwhm_pulse_duration/sqrt(2*log(2));  % Constant associated with the pulse width (seconds).
D_omega_tau = sqrt(-4 * log(chi) / T^2);

% Frequency vector steps
d_omega = 8 * D_omega_tau / number_of_samples; % Increment of the component frequencies (Hz).
% Calculate frequency differences
delta_omega = (-(number_of_samples / 2) * d_omega : d_omega : (number_of_samples / 2 - 1) * d_omega)';

% Calculate and normalize spectrum
ideal_spectrum = abs(exp(-delta_omega.^2/(4 * 2 * log(2) / fwhm_pulse_duration ^ 2))).^2; %exp(dw^2/4Gamma)
ideal_spectrum = ideal_spectrum / max(ideal_spectrum);
frequency_vector = delta_omega + carrier_angular_frequency;
ideal_spectrum((frequency_vector < cutoff_angular_frequency(1)) | (frequency_vector > cutoff_angular_frequency(2))) = [];
frequency_vector((frequency_vector < cutoff_angular_frequency(1)) | (frequency_vector > cutoff_angular_frequency(2))) = [];
        
fprintf('Generated pulse length = %e [fs]\n', 2 * pi * (2 * log(2) / pi) / fwhm(frequency_vector, ideal_spectrum) * 1e15);
wavelength_vector = (2 * pi * speed_of_light_vacuum ./ frequency_vector); 
ideal_spectrum_vector_in_wavelength = wrev(ideal_spectrum .* (2 * pi * speed_of_light_vacuum) ./ wavelength_vector.^2);
ideal_spectrum_vector_in_wavelength = ideal_spectrum_vector_in_wavelength / max(ideal_spectrum_vector_in_wavelength) + offset;
wavelength_vector = wrev(wavelength_vector);

fprintf('Start calculations for SNR = %i dB.\n', snr);

% Corrupt spectrum
corrupted_spectrum = awgn(ideal_spectrum_vector_in_wavelength, snr, 'measured');

% Process data
[processed_spectrum, intensity_time, ...
    angular_frequency_vector, time_vector, ...
    tau_from_omega_processed, ...
    tau_from_intensity_time, best_gauss_width, ...
    edge_window_value_full_process, window_width_full_process, ...
    error_val, recovered_SNR, ...
    recovered_carrier_angular_frequency, recovered_carrier_wavelength] = ...
    process_spectrum(wavelength_vector, corrupted_spectrum, ...
    speed_of_light_vacuum, initial_width_constant, supergauss_sigma, ...
    supergauss_pow, width_step,threshold, number_of_steps, SG_order, ...
    SG_window_size, remove_offset);
    
    error_theo = (tau_from_intensity_time - fwhm_pulse_duration)/fwhm_pulse_duration*100;
    
% Normalize
processed_spectrum = processed_spectrum/max(processed_spectrum);
processed_spectrum_wavelength = flipud(processed_spectrum * 2 * pi * speed_of_light_vacuum ./ (flipud(wavelength_vector).^2));
processed_spectrum_wavelength = processed_spectrum_wavelength / max(processed_spectrum_wavelength);

% Print results
fprintf('Pulse duration calculated from spectrum: %i [fs]\n', tau_from_omega_processed * 1e15);
fprintf('Pulse duration calculated from time intensity distribution: %i [fs]\n', tau_from_intensity_time * 1e15);