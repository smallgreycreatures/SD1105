x1 = 0; %[m]engine vertical position
x2 = 0; %[m]seismic mass vertical position
Fex = 0; %[N] The total vertical excitation force
m1 = 5000; %[kg] engine mass
m2 = 10000; %[kg] seismic mass
k1 = 2.5*10^6; %[N/m] spring rate of spring damper system between engine and seismic mass
k2 = 5*10^6; %[N/m] spring rate of the spring damper system between seismic mass and hull
c1 = 2000; %[Ns/m] coefficient of the spring-damper system betweeen engine and seismic mass
c2 = 5000; %[Ns/m] coefficient of the spring damper system between seismic mass and hull


file_id = fopen('kraft.bin', 'rb');
file = fread(file_id, 'float32');
frequency_data = file(1:2:length(file));
F0 = file(2:2:length(file));
index = 1;
max_value = 0;
disp('running')
disp(file(99,1))
fileID = fopen('test.txt','w');
fprintf(fileID,'%f\n',file);

%disp(frequency_data)

K = [k1, -k1;-k1, (k1+k2)]; %stiffness matrix (8a)
M = [m1, 0; 0, m2]; %mass matrix (8b)
C = [c1, -c1; -c1, (c1+c2)]; % damping matrix (8c)
length(F0)

%Eigen frequencies
eigen_freq1 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))-sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11a
eigen_freq2 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))+sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11b
    
disp('eigen')
disp(eigen_freq1)
disp(eigen_freq2)
hull_force_and_freq = [];
%max_hull_force_and_freq = [0,0];
max_hull_force_and_freq = calculate_max_hull_force_and_frequency(F0, frequency_data, K, C, M, eigen_freq1, eigen_freq2, k2, c2);
disp(max_hull_force_and_freq)

max_force_on_hull =  700;

disp('starting')
potential_candidates = [];
for m2 = 2000:1:3000
    M = [m1, 0; 0, m2]; %mass matrix (8b)


    %Variable Eigen frequencies
    eigen_freq1 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))-sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11a
    eigen_freq2 = (1/(2*pi))*sqrt((k1/(2*m1)+((k1+k2)/(2*m2)))+sqrt((k1/(2*m1))^2 +((k1+k2)/(2*m2))^2 + (k1^2 -(k1*k2))/(2*m1*m2))); %11b
    
    max_hull_force_and_freq = calculate_max_hull_force_and_frequency(F0, frequency_data, K, C, M, eigen_freq1, eigen_freq2, k2, c2);

    if max_hull_force_and_freq(1) < max_force_on_hull+1 && max_hull_force_and_freq(1) > max_force_on_hull-1
        disp('woop woop')
        disp(max_hull_force_and_freq)
        disp(m2)
        if round(real(max_hull_force_and_freq(1))) == max_force_on_hull
            disp('yippie')
        end
        potential_candidates = [potential_candidates; [max_hull_force_and_freq, m2]];
    end 
end
