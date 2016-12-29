function [ max_hull_force_and_freq ] = calculate_max_hull_force_and_frequency( F0, frequency_data, K, C, M, k2, c2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
index = 1;
hull_force_and_freq = [];
max_hull_force_and_freq = [0,0];
for F0_value = F0'
    
   f = frequency_data(index);

   ohm = 2*pi*f; %circular frequency
   %disp(f)
   matrix = (K + (i*ohm*C)-((ohm)^2*M)); %(Matrix for (8))
   
   %disp(matrix)
   
   F0_vector = [F0_value; 0];
   
   x = linsolve(matrix, F0_vector); %Equation (8)
   %disp(x)
  
    
    %disp(x(2))
    hull_force = (k2+ i*ohm*c2)*x(2);% (9)
    hull_force_and_freq = [hull_force_and_freq, [hull_force; f]];
    
    if hull_force > max_hull_force_and_freq(1)
        
        max_hull_force_and_freq = [hull_force, f];
        
    end
        
    index = index +1;
end

end
