function [r_new_i,VM_NEW_i,RFr_i,Ave_mach_i] =...
radial(i,N,M,mass_frac,r,Vm,gamma,R,T,delta_z)
%% Calculate Radial Changes
mass_frac(N,1:M) = 1;
mass_frac(1,1:M) = 0;

for j=1:N
y(j) = r(j,i);
x(j) = mass_frac(j,i);
if x(j) >= 1
x(j) = 1;
end
end

r_new_i(2:N-1,i) = pchip(x,y,mass_frac(2:N-1,1));


%Calculate Relaxation Factor for radial changes
tsum=0;
for j=1:N
Mach2_i= Vm(j,i) ^ 2 / (gamma * R * T(j,i));
tsum = tsum + Mach2_i;
end

Ave_mach_i = tsum / N;
clear tsum
RFr_i = 1/(1 + 5/24 * (1 - Ave_mach_i ^ 2) * ((r(N,i) - r(1,i)) / (delta_z(i))));

for j = 2:N-1;
    
r_new_i(j,i) = r(j,i) + 0.05 * (r_new_i(j,i) - r(j,i));  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

r_new_i(1,i) = r(1,i);
r_new_i(N,i) = r(N,i);

%% Calculate Vm Changes
for j = 1:N
y(j) = Vm(j,i);
x(j) = r(j,i);
end
VM_NEW_i(1:N,i) = pchip(x,y,r_new_i(1:N,i));