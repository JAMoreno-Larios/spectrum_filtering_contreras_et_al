function estimated_baseline = automatic_baseline_removal_shen_et_al(input_spectrum, threshold)
% Automatic baseline removal as described in 
% "Automatic baseline correction method for the open-path Fourier 
% transform infrared spectra by using simple iterative averaging" by
% Xianchun Shen et al.
% Vol. 26, No. 10 | 14 May 2018 | OPTICS EXPRESS A609 
% By José Agustín Moreno-Larios

    number_of_samples = length(input_spectrum);
    estimated_baseline = input_spectrum;
    previous_s_abs = 0;
    while(1)
        start_index = 1;
        end_index = number_of_samples - 2;
        while (start_index < floor(number_of_samples / 2))
            for index = start_index : end_index
                estimated_baseline(index + 1) = min(estimated_baseline(index + 1), (estimated_baseline(index) + estimated_baseline(index + 2) ) / 2);
            end
            start_index = start_index + 1;
            end_index = end_index - 1;
        end
        s_abs = sum(abs(input_spectrum - estimated_baseline));
        relative_error = abs(previous_s_abs - s_abs) / s_abs;
        if(relative_error < threshold || isnan(relative_error))
            break;
        end
        previous_s_abs = s_abs;
        input_spectrum = estimated_baseline;
    end
end