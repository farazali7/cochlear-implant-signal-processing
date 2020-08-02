clc
clear
close all

sampling_rate = 16000;
num_of_channels = 22; % Creates 22 channels with adjustment

% Generate bark scale domain frequency bands
b = hertz_to_bark_scale([100, 8000-0.0001]);
bark_vector = linspace(b(1), b(2), num_of_channels+1);
hertz_vector = bark_scale_to_hertz(bark_vector);

frequency_bands = cell(num_of_channels, 1);
for num = 1 : length(hertz_vector)-1
    if isempty(frequency_bands{1})
        overlap_factor = 1;
    else
        overlap_factor = 1;
    end 
    band = [hertz_vector(num)*overlap_factor hertz_vector(num+1)];
    frequency_bands{num} = band;
end

% Generate corresponding bandpass filters
[filters] = generate_bandpass_filters(frequency_bands, sampling_rate, 6);

% Read all the sound files in directory
sound_files = dir('original_sounds/*.wav');
for index=1:length(sound_files)
    audio_file = fullfile(sound_files(index).folder, sound_files(index).name);
    resampled_audio = read_and_resample(audio_file);
    
%     sound(resampled_audio, sampling_rate)
%     pause(size(resampled_audio, 1)/sampling_rate)
    
    [~, name, ext] = fileparts(audio_file);
    
    % Create channels for sound by filtering (use filtfilt to get minimum
    % phase digital filtering)
    filtered_channels = cell(length(filters), 1);
    for filter_num = 1 : length(filters)
        current_filter = filters{filter_num, 1};
        sos = current_filter{1, 1};
        gain = current_filter{1, 2};
        filtered_sound = filtfilt(sos, gain, resampled_audio);
        filtered_channels{filter_num, 1} = filtered_sound;
    end

%%
%     % Plot output signals of lowest and higher frequency channels
%     [num_samples, ~] = size(resampled_audio);
%     samples_vector = 1:num_samples;
%     
%     % Plot lowest frequency channel output
%     lowest_channel = filtered_channels{1, 1};
%     plot(samples_vector, lowest_channel)
%     xlabel("Sample Number");
%     ylabel("Amplitude");
%     title("Lowest Frequency Channel Output")
%     
%     freq_output_dir = "frequency_channel_outputs";
%     if ~exist(freq_output_dir, 'dir')
%         mkdir(freq_output_dir)
%     end
%     
%     low_freq_output = "frequency_channel_outputs/low_" + name + ".fig";
%     savefig(low_freq_output)
% 
%     % Plot highest frequency channel output
%     highest_channel = filtered_channels{16, 1};
%     plot(samples_vector, highest_channel)
%     xlabel("Sample Number");
%     ylabel("Amplitude");
%     title("Highest Frequency Channel Output")
%     
%     high_freq_output = "frequency_channel_outputs/high_" + name + ".fig";
%     savefig(high_freq_output)
%%
    % Get envelope of each frequency channel output
    envelopes = cell(length(filtered_channels), 1);
    rectified_signals = cell(length(filtered_channels), 1); % For later use
    for channel_num = 1:length(filtered_channels)
        % Rectify all signals
        rectified_signal = abs(filtered_channels{channel_num, 1});
        rectified_signals{channel_num, 1} = rectified_signal;
        
        % Low-pass filter the signals
        [zeroes, poles, gain] = cheby2(6, 20, 400/(sampling_rate/2));
        [sos, gain] = zp2sos(zeroes, poles, gain);
        envelope = filtfilt(sos, gain, rectified_signal);
        envelopes{channel_num, 1} = envelope;
    end

%%
%     % Plot envelopes of lowest and higher frequency channel outputs
% 
%     % Plot lowest frequency channel envelope
%     lowest_channel_envelope = envelopes{1, 1};
%     plot(samples_vector, lowest_channel_envelope)
%     xlabel("Sample Number");
%     ylabel("Amplitude");
%     title("Lowest Frequency Channel Envelope")
%     
%     envelope_dir = "frequency_channel_envelopes";
%     if ~exist(envelope_dir, 'dir')
%         mkdir(envelope_dir)
%     end
%     
%     low_freq_envelope = "frequency_channel_envelopes/low_" + name + ".fig";
%     savefig(low_freq_envelope)
% 
%     % Plot highest frequency channel envelope
%     highest_channel_envelope = envelopes{16, 1};
%     plot(samples_vector, highest_channel_envelope)
%     xlabel("Sample Number");
%     ylabel("Amplitude");
%     title("Highest Frequency Channel Envelope")
%     
%     high_freq_envelope = "frequency_channel_envelopes/high_" + name + ".fig";
%     savefig(high_freq_envelope)
    
    %%
    amplitude_modulated_cosine_signals = cell(length(filtered_channels), 1);
    % Generate cosine signal with each channels' central frequency
    for channel_num = 1:length(frequency_bands)
        band_low = frequency_bands{channel_num}(1);
        band_high = frequency_bands{channel_num}(2);
        
        ratio = band_high / band_low;
        if ratio >= 1.1
            center_frequency = sqrt(band_high * band_low);
        else
            center_frequency = (band_high + band_low)/2;
        end
        
        dt = 1/sampling_rate;
        audio_length = size(resampled_audio, 1);
        duration = (0:audio_length-1)/sampling_rate;
        tt = 0:dt:duration-1;
        
        cos_signal = cos(2*pi*center_frequency*duration);
        
        cos_signal = cos_signal.';
        
        % Amplitude modulated cosine signal
        rectified_signal = envelopes{channel_num, 1};
        amplitude_modulated_signal = cos_signal .* rectified_signal;
        
        amplitude_modulated_cosine_signals{channel_num, 1} = amplitude_modulated_signal; 
    end
    
    % Add all signals together
    output_signal = zeros(size(resampled_audio));
    for signal_num = 1:size(amplitude_modulated_cosine_signals, 1)
        current_signal = amplitude_modulated_cosine_signals{signal_num, 1};       
        output_signal = output_signal + current_signal;
    end
    
    % Normalize output signal
    output_signal_normalized = output_signal / max(abs(output_signal));
     
    % Play output signal
%     sound(output_signal_normalized, sampling_rate)
%     pause(audio_length/sampling_rate)
    
%     % Write sound to new file
%     output_sound_files_dir = "prototype3_output_sound_files";
%     if ~exist(output_sound_files_dir, 'dir')
%         mkdir(output_sound_files_dir)
%     end
%         
%     filename = "prototype3_output_sound_files/output_" + name + ".wav";
%     audiowrite(filename, output_signal_normalized, sampling_rate);

end

