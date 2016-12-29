%Written by Frans Nordén 
%Last changed 2016-12-29


function [ max_hull_force_and_freq ] = calculate_max_hull_force_and_frequency( F0, frequency_data, K, C, M, k2, c2)
%calculate_max_hull_force_and_frequency
%   loops through a list of amplitude and frequencies
%   and decides based of formulas what the greatest force are.
%   The greatest force are returned togeteher with its matching frequency
index = 1;
max_hull_force_and_freq = [0,0];
for F0_value = F0' %transposed for easier access
    
   f = frequency_data(index);

   ohm = 2*pi*f; %circular frequency

   matrix = (K + (i*ohm*C)-((ohm)^2*M)); %Matrix for (8)
   
   F0_vector = [F0_value; 0]; %(8e)
   
   x = linsolve(matrix, F0_vector); %Solving Equation 8

    hull_force = (k2+ i*ohm*c2)*x(2);% (9)
    
    if hull_force > max_hull_force_and_freq(1)
        
        max_hull_force_and_freq = [hull_force, f];
        
    end
        
    index = index +1;
end

end
