%Written by Frans Nordén 
%Last changed 2016-12-29
clear

x1 = 0; %[m]engine vertical position
x2 = 0; %[m]seismic mass vertical position
Fex = 0; %[N] The total vertical excitation force
m1 = 5000; %[kg] engine mass
m2 = 10000; %[kg] seismic mass
k1 = 2.5*10^6; %[N/m] spring rate of spring damper system between engine and seismic mass
k2 = 5*10^6; %[N/m] spring rate of the spring damper system between seismic mass and hull
c1 = 2000; %[Ns/m] coefficient of the spring-damper system betweeen engine and seismic mass
c2 = 5000; %[Ns/m] coefficient of the spring damper system between seismic mass and hull

min_mass = 2000; %[kg] start value for searching for suitable seismic mass
max_mass = 3000; %[kg] stop value for searching for suitable seismic mass
max_force_on_hull = 700; %[N] max force as said in requirments
step = 1; %how many kilos to step at a time in find_suitable_mass

file_id = fopen('kraft.bin', 'rb');
file = fread(file_id, 'float32');
frequency_data = file(1:2:length(file));
F0 = file(2:2:length(file));

%disp(frequency_data)

K = [k1, -k1;-k1, (k1+k2)]; %stiffness matrix (8a)
M = [m1, 0; 0, m2]; %mass matrix (8b)
C = [c1, -c1; -c1, (c1+c2)]; % damping matrix (8c)

%Eigen frequencies (11a and 11b)
eigen_freq1 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))-sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11a
eigen_freq2 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))+sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11b
    
max_hull_force_and_freq = calculate_max_hull_force_and_frequency(F0, frequency_data, K, C, M, k2, c2);
if max_hull_force_and_freq(2) <= eigen_freq1 || max_hull_force_and_freq(2) <= eigen_freq2
    disp(['Max hull force is ',num2str(max_hull_force_and_freq(1)), 'N and frequency=', num2str(max_hull_force_and_freq(2)), 'Hz while eigen frequencies are ', num2str(eigen_freq1), 'Hz and ', num2str(eigen_freq2), 'Hz'])
    disp('This does not fullfill the requirments')
    
    suitable_mass = find_suitable_mass( min_mass, max_mass, step, m1, F0, frequency_data, K, C, k2, c2, max_force_on_hull);
   
    disp(['To fullfill the requirments the seismic mass should be ', num2str(suitable_mass), 'kg'])
else
    disp(['Max hull force is ',num2str(max_hull_force_and_freq(1)), 'N and frequency=', num2str(max_hull_force_and_freq(2)), 'Hz while eigen frequencies are ', num2str(eigen_freq1), 'Hz and ', num2str(eigen_freq2), 'Hz'])
    disp('This does fullfill the requirments')
end


