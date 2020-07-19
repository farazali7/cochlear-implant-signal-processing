% [10] H. Traunmüller, “Analytical expressions for the tonotopic sensory scale,” The Journal of the Acoustical Society of America, vol. 88, no. 1, pp. 97–100, 1990. 

function [barks] = hertz_to_bark_scale(hertz)
    barks = [];
    for num = 1 : length(hertz)
        curr_hertz = hertz(num);
        numerator = 26.81 * curr_hertz;
        denominator = 1960 + curr_hertz;
        bark = (numerator / denominator) - 0.53;
        % Value correction
        if bark < 2
            bark = bark + (0.15 * (2 - bark));
        elseif bark > 20.1
            bark = bark + (0.22 *(bark - 20.1));
        end
        
        barks(num) = bark;
    end
       
end

