function audio_resampled = read_and_resample(x)
    % Get file name 
    [~, name, ext] = fileparts(x);
    
    % Read and find sample rate of audio
    [audio,sample_rate] = audioread(x);
    
    % Stereo conversion to mono
    [num_samples,channels] = size(audio);
    
    if channels == 2
        audio_mono = sum(audio, 2) / size(audio, 2);
    else 
        audio_mono = audio;
    end
    
%     % Play sound
%     sound(audio_mono, sample_rate);
%     pause(num_samples/sample_rate)
%     
%     % Write sound to new file
%     mono_sound_files_dir = "mono_sound_files";
%     if ~exist(mono_sound_files_dir, 'dir')
%         mkdir(mono_sound_files_dir)
%     end
%         
%     filename = "mono_sound_files/mono_" + name + ".wav";
%     audiowrite(filename, audio_mono, sample_rate);
    
%     % Plot sound waveform as function of sample number
%     orig_samples = (1:num_samples);
%     
%     plot(orig_samples, audio_mono);
%     xlabel("Sample Number");
%     ylabel("Amplitude");
%     
%     waveforms_dir = "waveforms";
%     if ~exist(waveforms_dir, 'dir')
%         mkdir(waveforms_dir)
%     end
%     
%     old_waveform = "waveforms/original_" + name + ".fig";
%     savefig(old_waveform)
        
    % Resample input signal to 16kHz
    if (sample_rate >= 16000)
        audio_resampled = resample(audio_mono, 16000, sample_rate);
    else
        error("Sample rate too low for audio file: " + name + ext)
    end
    
end