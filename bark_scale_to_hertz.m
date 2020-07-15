function [hertz] = bark_scale_to_hertz(barks)
    hertz = [];
    for num = 1 : length(barks)
        curr_bark = barks(num);
        % Value correction
        if curr_bark < 2
            curr_bark = (curr_bark - 0.3)/0.85;
        elseif curr_bark > 20.1
            curr_bark = (curr_bark + 4.422)/1.22;
        end
        
        numerator = curr_bark + 0.53;
        denominator = 26.28 - curr_bark;
        hz_value = (numerator/denominator) * 1960;
        hertz(num) = hz_value;
    end

end

