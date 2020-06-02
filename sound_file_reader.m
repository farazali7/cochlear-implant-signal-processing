clear all
close all


function output = read_and_resample(x)
    % Read and find sample rate of audio
    [audio,sample_rate] = audioread(x);
    
    % Stereo conversion to mono
    [num_samples,channels] = size(audio);
    
    if n == 2
        audio_mono = sum(audio, 2) / size(audio, 2);
    else 
        audio_mono = audio;
    end
    
    % Play sound
    sound(audio_mono, sample_rate);
    
    % Write sound to new file - TODO: location/names
    filename = "new_sound.wav";
    audiowrite(filename, audio_mono, sample_rate);
    
    % Plot sound waveform as function of sample number
    time = (1:num_samples)/sample_rate;
    plot(time, audio_mono);
    xlabel("Time");
    ylabel("Amplitude");
    
    % Resample input signal to 16kHz - TODO: lower than 16kHz case
    audio_resampled = resample(audio_mono, 16000, sample_rate);
    
    % Generate signal with cosine.....? confused idk about this one
    
        
end