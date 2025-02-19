function bandpass_filters = generate_bandpass_filters(frequency_bands, sampling_rate, filter_order)
    bandpass_filters = cell(length(frequency_bands), 1);
    for num = 1 : length(frequency_bands)
        band = frequency_bands{num};
        normalized_band = band/(sampling_rate/2);
        
        % Chebyshev Type II Filter
        [zeros, poles, gain] = cheby2(filter_order/2, 20, normalized_band);
        [sos, gain] = zp2sos(zeros, poles, gain);
        filter_cell = {sos, gain};
        bandpass_filters{num, 1} = filter_cell;
    end
end