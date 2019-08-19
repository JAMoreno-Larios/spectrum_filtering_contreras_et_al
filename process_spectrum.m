function [processed_spectrum, intensity_time, ...
    angular_frequency_vector, time_vector, ...
    tau_from_omega_processed, tau_from_intensity_time, best_width_constant, ...
    edge_window_value_full_process, window_width_full_process, ...
    error_val, SNR, carrier_angular_frequency, carrier_wavelength] = ...
    process_spectrum(wavelength_vector, spectrum_vector_in_wavelength, ...
    speed_of_light_vacuum, initial_width_constant, supergauss_sigma, ...
    supergauss_pow, width_step,threshold, number_of_steps, ...
    SG_order, SG_window_size, remove_offset)

    % Get number of samples
    number_of_samples = length(spectrum_vector_in_wavelength);
    fft_samples = 4 * number_of_samples;
    final_number_of_samples = 2 * number_of_samples;
    
    % Define frequency and time vector
    angular_frequency_vector = wrev(2 * pi * ...
        speed_of_light_vacuum ./ wavelength_vector);
    freq_step = mean(diff(angular_frequency_vector));

    % Prepare for fft
    time_step = 2 * pi / (fft_samples * freq_step);
    time_vector = (-fft_samples / 2 : 1 : ...
        fft_samples / 2 - 1)' * time_step;
    time_vector((fft_samples + final_number_of_samples)/2 +1 : end) = [];
    time_vector(1 : (fft_samples - final_number_of_samples)/2) = [];
    
    % Apply Savitzky - Golay filter
    order = SG_order;
    framelen = SG_window_size;
    spectrum_vector_in_wavelength_smoothed = ...
        sgolayfilt(spectrum_vector_in_wavelength, order, framelen);
    
    % Get noise
    noise = spectrum_vector_in_wavelength - ...
        spectrum_vector_in_wavelength_smoothed;

    % Get power
    P_signal = (1/number_of_samples) * ...
        sum(spectrum_vector_in_wavelength_smoothed.^2);
    P_noise = (1/number_of_samples) * sum(noise.^2);

    % Get SNR
    SNR = 10*log10(P_signal / P_noise);
    
    % Remove offset
    if(remove_offset)
        threshold_baseline = 0.01;
        spectrum_wavelengtht_baseline = ...
            automatic_baseline_removal_shen_et_al...
            (spectrum_vector_in_wavelength_smoothed,threshold_baseline);
        spectrum_vector_in_wavelength_processed = ...
            spectrum_vector_in_wavelength_smoothed ...
            - spectrum_wavelengtht_baseline;
    else
        spectrum_vector_in_wavelength_processed = ...
            spectrum_vector_in_wavelength_smoothed;
    end
    
    spectrum_vector_in_wavelength_processed = ...
        spectrum_vector_in_wavelength_processed / ...
        max(spectrum_vector_in_wavelength_processed);

    % Convert S(lambda) to S(omega) - apply necessary flip and renormalize
   spectrum_frequency = wrev(spectrum_vector_in_wavelength_processed .* ...
        wavelength_vector .^ 2 / ...
        (2 * pi * speed_of_light_vacuum));
    spectrum_frequency = spectrum_frequency / max(spectrum_frequency);
    
    % Get carrier wavelength via centroid
    sum_of_spectrum_elements = sum(spectrum_frequency);
    carrier_angular_frequency = sum(angular_frequency_vector .* ...
        spectrum_frequency / sum_of_spectrum_elements);
    carrier_wavelength = 2 * pi * ...
        speed_of_light_vacuum / carrier_angular_frequency;

    % Get delta_omega
    delta_omega = fwhm(angular_frequency_vector, spectrum_frequency);
    % Estimate delta_tau
    tau_from_omega_after_sg = 2 * pi * (2 * log(2) / pi) / delta_omega;
    
    % Initialize best window search algorithm
    width_constant = initial_width_constant;
    % Get width
    frequency_fwhm = delta_omega;
           
    index = 1;
    while(index <= number_of_steps)

    % Apply supergauss
        window_amplitude = ...
            supergauss_window(angular_frequency_vector, ...
            width_constant * frequency_fwhm, carrier_angular_frequency, ...
            supergauss_sigma, supergauss_pow);
        processed_spectrum = ...
            spectrum_frequency .* ...
            window_amplitude;

        % Build electric field
        initial_phase = 0;
        electric_field_frequency = sqrt(processed_spectrum) .* ...
            exp(1i * initial_phase);
        
        % Pad with zeros to improve time resolution
        electric_field_frequency = padarray(electric_field_frequency, ...
            (fft_samples - number_of_samples) / 2);

        % Get electric field in time
        electric_field_time = ifftshift(ifft(electric_field_frequency));

        % Get intensity in time and normalize
        intensity_time = abs(electric_field_time) .^ 2;
        intensity_time = intensity_time / max(intensity_time);
        
        % Remove extra components
        intensity_time((fft_samples + final_number_of_samples)/2 +1 : end) = [];
        intensity_time(1 : (fft_samples - final_number_of_samples)/2) = [];

        % Get FWHM from temporal intensity
        tau_from_intensity_time = fwhm(time_vector, intensity_time);

        % Calculate error
        error_val = (tau_from_intensity_time - ...
            tau_from_omega_after_sg)/tau_from_omega_after_sg*100;
        best_width_constant = width_constant;
        window_width_full_process = width_constant * frequency_fwhm;
        edge_window_value_pos = find(window_amplitude, 1, 'last');
        if(isempty(edge_window_value_pos))
            warning('Could not measure FWHM at intensity')
            edge_window_value_full_process = NaN;
            break;
        else
            edge_window_value_full_process =  ...
                window_amplitude(edge_window_value_pos);
        end
        if(abs(error_val) <= threshold)
            break;
        elseif(error_val < 0)
            width_constant = width_constant - width_step;
        elseif(error_val > 0)
            width_constant = width_constant + width_step;
        end
        index = index + 1;
        if(index == number_of_steps)
            warning('Could not get an optimal width value using the stated parameters')
        end
    end
    
    % Measure carrier again
    sum_of_spectrum_elements = sum(processed_spectrum);
    carrier_angular_frequency = sum(angular_frequency_vector .* ...
        processed_spectrum / sum_of_spectrum_elements);
    carrier_wavelength = 2 * pi * ...
        speed_of_light_vacuum / carrier_angular_frequency;
    
    % Get delta_omega processed
    delta_omega_processed = fwhm(angular_frequency_vector, ...
        processed_spectrum);
    % Estimate delta_tau
    tau_from_omega_processed = 2 * pi * ...
        (2 * log(2) / pi) / delta_omega_processed;
end