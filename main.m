clc
clear
close all

sampling_rate = 16000;
num_of_channels = 15; % Creates 16 channels with adjustment

% Generate bark scale domain frequency bands
b = hertz_to_bark_scale([100, 6700]);
bark_vector = linspace(b(1), b(2), num_of_channels+1);
hertz_vector = bark_scale_to_hertz(bark_vector);

frequency_bands = cell(num_of_channels+1, 1);
for num = 1 : length(hertz_vector)-1 
    band = [hertz_vector(num) hertz_vector(num+1)];
    frequency_bands{num} = band;
end

% Add in last channel
last_band_low = frequency_bands{num_of_channels}(2);
last_band_high = double(8000) - (0.0001); % Adjust for normalized Nyquist Cutoff Freq
last_band = [last_band_low last_band_high];
frequency_bands{num_of_channels+1} = last_band;

% Generate corresponding bandpass filters
[filters] = generate_bandpass_filters(frequency_bands, sampling_rate, 6);

% Read all the sound files in directory
sound_files = dir('original_sounds/*.wav');
for index=1:length(sound_files)
    audio_file = fullfile(sound_files(index).folder, sound_files(index).name);
    resampled_audio = read_and_resample(audio_file);
    
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
        [zeroes, poles, gain] = butter(6, 400/(sampling_rate/2));
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
         rectified_signal = rectified_signals{channel_num, 1};
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
    sound(output_signal_normalized, sampling_rate)
    pause(audio_length/sampling_rate)
    
%     % Write sound to new file
%     output_sound_files_dir = "output_sound_files";
%     if ~exist(output_sound_files_dir, 'dir')
%         mkdir(output_sound_files_dir)
%     end
%         
%     filename = "output_sound_files/output_" + name + ".wav";
%     audiowrite(filename, output_signal_normalized, sampling_rate);

end

%     % Generate cosine signal with 1kHz oscillations and same time duration
%     % as original audio
%     fs = 16000; % Sampling frequency (samples per second) 
%     F = 1000; % Cosine wave frequency (hertz) 
%     
%     dt = 1/fs;
%     
%     % Original time of audio in seconds
%     time_duration = size(audio_resampled,1)/fs;
%     tt = 0:dt:time_duration;
%     
%     cos_signal_sound = cos(2*pi*F*tt);
%     
%     % Play sound generated by cosine signal
%     sound(cos_signal_sound, fs)
%     pause(time_duration)
%     
%     % Plot two cycles of cosine sound signal waveform w.r.t time
%     T = 2/F;
%     
%     tt_plot = 0:dt:T;
%     cos_signal_plot = cos(2*pi*F*tt_plot);
%         
%     plot(tt_plot,cos_signal_plot);
%     xlabel("Time")
%     ylabel("Amplitude");
%     
%     cosine_waveform = "waveforms/cosine_" + name + ".fig";
%     savefig(cosine_waveform)

