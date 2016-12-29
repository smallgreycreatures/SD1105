%Written by Frans Nordén 
%Last changed 2016-12-29


function [ calculated_mass ] = find_suitable_mass( min_mass, max_mass, step, m1, F0, frequency_data, K, C, k2, c2, max_force_on_hull)
%find_suitable_mass
%   Find a suitable mass for the seismic mass through testing what 
%   values that might be suitable, i.e every kilogram between min_mass and
%   max_mass are tested and inserted into
%   calculate_max_hull_force_and_frequency. If a suitable mass is found it
%   is returned. If no mass is found it tries with better precision.

calculated_mass = 0;
step_changer = 1; %Factor to scale down the step taken 10^-(step_changer)
for i = 1:8 %seems unreasonable to go under 10^-8
    for m2 = min_mass:step:max_mass
        M = [m1, 0; 0, m2]; %mass matrix (8b)

        max_hull_force_and_freq = calculate_max_hull_force_and_frequency(F0, frequency_data, K, C, M, k2, c2);

        if max_hull_force_and_freq(1) < max_force_on_hull+1 && max_hull_force_and_freq(1) > max_force_on_hull-1

            if round(real(max_hull_force_and_freq(1))) == max_force_on_hull %rounds to closest integer as in req
                calculated_mass = m2;
            end
        
        end 
    end
    if calculated_mass ~= 0
       return
    else
        step = step*10^(-step_changer); %changes step to step*10^-i
    end
end
end

