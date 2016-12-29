function [ calculated_mass ] = find_suitable_mass( min_mass, max_mass, m1, F0, frequency_data, K, C, k2, c2, max_force_on_hull)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

calculated_mass = 0;
for m2 = min_mass:1:max_mass
    M = [m1, 0; 0, m2]; %mass matrix (8b)

    max_hull_force_and_freq = calculate_max_hull_force_and_frequency(F0, frequency_data, K, C, M, k2, c2);

    if max_hull_force_and_freq(1) < max_force_on_hull+1 && max_hull_force_and_freq(1) > max_force_on_hull-1

        if round(real(max_hull_force_and_freq(1))) == max_force_on_hull
            calculated_mass = m2;
        end
        
    end 
end

end

