function [SM_new_i,VM_NEW_i] =...
radial(i,N,M,mass_frac,SM,Vm,gamma,R,T,delta_z,r)
%% Calculate Radial Changes
mass_frac(N,1:M) = 1;
mass_frac(1,1:M) = 0;

for j=1:N
y(j) = SM(j,i);
x(j) = mass_frac(j,i);
if x(j) >= 1
x(j) = 1;
end
end

SM_new_i(2:N-1,i) = pchip(x,y,mass_frac(2:N-1,1));


for j = 2:N-1;
    
SM_new_i(j,i) = SM(j,i) + 0.5 * (SM_new_i(j,i) - SM(j,i));  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

SM_new_i(1,i) = SM(1,i);
SM_new_i(N,i) = SM(N,i);

%% Calculate Vm Changes
for j = 1:N
y(j) = Vm(j,i);
x(j) = SM(j,i);
end
VM_NEW_i(1:N,i) = pchip(x,y,SM_new_i(1:N,i));