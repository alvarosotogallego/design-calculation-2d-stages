%% Main
clear
clc
close all
tic
%Inputs

global r z bladeTest TIPO
global epsilon psi rel_f
global gastoM presion0 temperatura0 omegaMatriz espesorMatriz omega
global PerdidaCT PerdidaCH PerdidaOMEGA 
global PB12 PB22 PB32 PB42 PB12P PB22P PB32P PB42P Ticks02 TicksBeta2 


%Velocity Tolerance
Vel_tol=4;
%Mass Flow Rate Tolerance
mass_flow_tol = 1e-1;
%Velocity Relaxation Factor
RFv = 0.1;
%Gas Constant
R = 287;
%Specific Heat (Constant Pressure)
Cp = 1004.5;
%Specific Heat Ratio
gamma = 1.4;

V_tREL = zeros([50 100]); 
V_tREL1 = zeros([50 100]);
V_ang = zeros([50 100]);

%Number of Streamlines
N = 50;
%Compute the number of Axial Locations
M = 100;

z_h(1:M)=z(1,:);
z_c(1:M)=z(end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constant Input Stagnation Pressure
p_0_In = presion0;
%Constant Input Stagnation Temperature
T_0_In = temperatura0;
%Radial Velocity at Inlet
Vr_In=0;

%Axial Velocity at Inlet (Estimación)
Vz_In = 100;
test =0;
%Blade Angles (Unbladed Region - Blade angles are zero)

%Annulus Loss Coefficients
omega_m=PerdidaOMEGA;
ct=PerdidaCT;
ch=PerdidaCH;

%Total Velocity at Inlet (Estimación)
V_In = sqrt( Vz_In ^ 2 + Vr_In ^ 2 );
%Constat Static Temperature at Inlet (Estimación)
T_In(:,1) = T_0_In(:,1) - V_In(:,1).^ 2 / (2 * Cp);
%Entropy at Inlet (Estimación)
s_In = Cp * log(T_0_In) - R * log(p_0_In);
%Density at Inlet (Estimación)
rho_In = (T_In ^ ((Cp - R) / R) * exp(-s_In / R)) / R;

%Axial Velocity at Inlet (Estimación)
Vz_In = gastoM/(rho_In * 2 * pi * (z_c(1)-z_h(1)) * r(1,1) )
%Total Velocity at Inlet (Estimación)
V_In = sqrt( Vz_In ^ 2 + Vr_In ^ 2 );
%Constat Static Temperature at Inlet (Estimación)
T_In(:,1) = T_0_In(:,1) - V_In(:,1).^ 2 / (2 * Cp);
%Entropy at Inlet (Estimación)
s_In = Cp * log(T_0_In) - R * log(p_0_In);
%Density at Inlet (Estimación)
rho_In = (T_In ^ ((Cp - R) / R) * exp(-s_In / R)) / R;

%Axial Velocity at Inlet (Real)
Vz_In = gastoM/(rho_In * 2 * pi * (z_c(1)-z_h(1)) * r(1,1) )
%Total Velocity at Inlet (Real)
V_In = sqrt( Vz_In ^ 2 + Vr_In ^ 2 );
%Constat Static Temperature at Inlet (Real)
T_In(:,1) = T_0_In(:,1) - V_In(:,1).^ 2 / (2 * Cp);
%Entropy at Inlet (Real)
s_In = Cp * log(T_0_In) - R * log(p_0_In);
%Density at Inlet (Real)
rho_In = (T_In ^ ((Cp - R) / R) * exp(-s_In / R)) / R;

%Mass flow at Itlet
mdot_tot = gastoM;
%Mass flow per steamline
mdot_per = mdot_tot / (N - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:M
delta_z(i) = z(i) - z(i - 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Middle Streamline Index Number
SL_Mid = round( N / 2);

%Tangential Velocity Profile
for i=1:M
for j=1:N
            V_t(j,i)=0;
end
end


Vz_In(1:N) = Vz_In;
Vz(1:N,1) = Vz_In;

c1u=0;

%Use simplified radial equilibrium to find Vz across inlet
%Below Midline
for j = (SL_Mid - 1):-1:1
Vz_In(j) = sqrt(Vz_In(j+1) ^ 2 - 1/2 * (r(j,1) - r(j+1,1)) * (1 / (r(j,1) ^ 2) * deriv(j,(r(:,1).^ 2.* V_t(:,1).^ 2),r(:,1),N,1) + 1 / (r(j+1,1) ^ 2) * deriv(j + 1,(r(:,1).^ 2.* V_t(:,1).^ 2),r(:,1),N,1)));
end
%Above Midline
for j = SL_Mid + 1:N

Vz_In(j) = sqrt(Vz_In(j-1) ^ 2 - 1/2 * (r(j,1) - r(j-1,1)) * (1 / (r(j,1) ^ 2) * deriv(j,(r(:,1).^ 2.* V_t(:,1).^ 2),r(:,1),N,1) + 1 / (r(j-1,1) ^ 2) * deriv(j-1,(r(:,1).^ 2.* V_t(:,1).^ 2),r(:,1),N,1)));
end

%Proportion of the mass flow rate in each stream tube
mass_frac(1) = 0;
for j = 2:N
mass_frac(j) = mdot_per / mdot_tot + mass_frac(j - 1);
end

%Calculate the Axial Velocity Based on Mass Flow rates and Initial Streamlines
for i = 2:M
Vz(1:N,i) = mdot_tot / (rho_In * pi *    (  (r(N,i))^2 - (r(1,i))^2  +  r(1,i)*2*(z(N,i)-z(1,i))  )  );
end

%Extend initial Pressures across annulus for first estimate
for j = 1:N
for i = 1:M
p_0(j,i) = p_0_In;
end
end

Vm = Vz;

%Update Properties
H_0_In(1:N) = Cp * T_0_In;
H_In = H_0_In - .5*V_In.^2;
T_0(1:N,1:M) = T_0_In;
V_tot = sqrt( Vm.^2 + V_t.^2 );
T = T_0 - (V_tot.^2) / (2 * Cp);
H_0(1:N,1:M) = H_0_In(1);
H = Cp * T_0-(V_tot.^2) / 2;
s = Cp * log(T_0) - R * log(p_0);
p = exp((Cp./ R.* log(T))-s./ R);
rho = p./ ( R.*T );

for k=1:N
mass_f(k,1) = mass_frac(k);
end

H_In=H_In(1);

mdot_per = mdot_per * 1.01;

for j=2:N-1
r(j,end) = sqrt( mdot_per / (rho(j,end) * Vm(j,end) * pi) + r(j - 1,end) ^ 2);
end

for i=1:M-1
for j=2:N-1 
r(j,i) = (r(j,end) - r(1,end)) / (r(N,end) - r(1,end)) * (r(N,i)-r(1,i)) + r(1,i);
end
end

mdot_perR2 = zeros([50 100]);
mdot_TOT2 = 0;

for i=1:M
for j=2:N
    
mdot_perR2(j,i) = (r(j,i)^2 - r(j-1,i)^2 + r(1,i)*2*(z(j,i)-z(j-1,i)) )* rho(j,i) * Vz(j,i) * pi;

if i== 100
    
    mdot_TOT2 = mdot_TOT2 + mdot_perR2(j,i);
    
end
end
end

mdot_perR2 

for i = 1:M
for j=1:N
    if i == 1
        V_t(j,i) = omega*2*3.1416/60*r(j,i);
    elseif bladeTest(j,i) == 2
        V_t(j,i)=tan(rel_f(j,i))*Vm(j,i);
    else      
        V_t(j,i)= V_t(j,i-1)*r(j,i-1)/r(j,i); 
    end
end
end

%clear variables that have been restructured into other variables
clear T_0_In T_In H_0_In p_In p_0_In
clear V_t_In rho_In s_In r_ft r_h r_c V_In Vr_In
clear mass_frac
mass_frac = mass_f;
clear mass_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% BEGINNING OF OUTER ITERATION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Mass Iteration Counter Equal to Zero
mass_iter(1:M) = 0;
%Set Velocity Iteration Counter Equal to Zero
V_iter(1:M) = 0;
%Set Radial Change Iteration Counter Equal to Zero
r_iter(1:M) = 0;
%Set Velocity and Radial Position Flags
V_flag = 1;
R_flag = 1;

% while R_flag == 1

while V_flag == 1
%% Curvature and Flow Angle
for i = 1:M
[C_i,phi_i] = curve(i,r,z,N,M);
% C(1:N,i) = C_i(1:N,i);
C(1:N,i) = 0;
phi(1:N,i) = phi_i(1:N,i);
end

Vz = Vm.*cos(phi);
Vr = Vm.*sin(phi);

% pcolor(z,r,C); 
% colorbar;
% colormap('jet');
% title('C');

if mass_iter(i) > 0

Vm(1:N,1) = Old_Vm(1:N,1);
end
clear phi_i


%% Gradient
for i = 1:M
for j = 1:N %units of ft
if i == M
delta_m(j,M) = delta_m(j,M-1);
rV_t2(j,i) = (r(j,i) * V_t(j,i)) ^ 2;
else
delta_m(j,i) = sqrt(((r(j,i+1) - r(j,i))) ^ 2 + ...
(z(i+1) - z(i)) ^ 2);

rV_t2(j,i) = (r(j,i) * V_t(j,i)) ^ 2;
end
end
end

for i = 1:M
[rV_t2_i,s_i,V_t_i,T_0_i,H_i,T_i,p_i,p_0_i,V_tot_i,rho_i,dsdq_i,...
drvtdr_i,dHodq_i,dvmdm_i,delta_s_i,delta_po_i] = ...
gradient(i,r,T_0,N,Vm,rho,V_t,s,H_0,Cp,R,gamma,T,SL_Mid,p_0,V_tot,...
p,delta_m,rV_t2,M,omega_m,ct,ch,H,test);
rV_t2(1:N,i) = rV_t2_i(1:N,i);
s(1:N,i) = s_i(1:N,i);
V_t(1:N,i) = V_t_i(1:N,i);
T_0(1:N,i) = T_0_i(1:N,i);
H(1:N,i) = H_i(1:N,i);
T(1:N,i) = T_i(1:N,i);
p(1:N,i) = p_i(1:N,i);
p_0(1:N,i) = p_0_i(1:N,i);
V_tot(1:N,i) = V_tot_i(1:N,i);
rho(1:N,i) = rho_i(1:N,i);
dsdq(1:N,i) = dsdq_i(1:N,i);
drvtdr(1:N,i) = drvtdr_i(1:N,i);
dHodq(1:N,i) = dHodq_i(1:N,i);
dvmdm(1:N,i) = dvmdm_i(1:N,i);
delta_s(1:N,i) = delta_s_i(1:N,i);
delta_po(1:N,i) = delta_po_i(1:N,i);
end

clear rV_t2_i s_i H_0_i V_t_i T_0_i H_i T_i p_i p_0_i V_tot_i rho_i
clear k_ht_i dsdq_i drvtdr_i dHodq_i dvmdm_i delta_s_i delta_po_i
%% Velocity
test =1;
Old_Vm = Vm;
for i = 2:M
mass_frac(1:N,i) = mass_frac(1:N,1);
end

%Calculate blade forces (no hacen nada)
for i = 2:M
for j = 1:N
num(j,i) = r(j,i) * Vm(j,i) * tan(rel_f(j,i)) - r(j,i-1) * Vm(j,i-1) * tan(rel_f(j,i-1));
delta_m(j,i) = sqrt( ((r(j,i) - r(j,i-1))) ^ 2 + (z(j,i) - z(j,i-1)) ^ 2);
d1dm(j,i) = num(j,i) / delta_m(j,i);
end
end
for i = 1:M
for j = 1:N
FN(j,i) = Vm(j,i) / r(j,i) * (tan(rel_f(j,i)) * tan(phi(j,i) + psi(j,i)) + tan(epsilon(j,i)) / cos(phi(j,i) + psi(j,i))) * d1dm(j,i);
FM(j,i) = Vm(j,i) * tan(rel_f(j,i)) / r(j,i) * d1dm(j,i);
end
end

%s

for i = 2:M
  
    
[Vm_i,mass_frac_i,mass_iter_i,MASSFLOW,A_i,B_i,V_t,T_0,T,p,p_0,s,H,H_0,V_tot,V_tREL1,V_ang,V_tREL,rho] = ...
velocity_TEST(i,SL_Mid,Vm,phi,C,N,dvmdm,dHodq,dsdq,drvtdr,r,RFv,rho,...
mdot_tot,mass_flow_tol,mass_iter,FM,FN,T_0,V_t,bladeTest,rel_f,omegaMatriz,Cp,gamma,T,p,p_0,s,H,H_0,R,V_tot,V_tREL1,V_tREL,V_ang,TIPO,z);

Vm_NEW(1:N,i) = Vm_i(1:N,i);
mass_frac(1:N,i) = mass_frac_i(1:N,i);
mass_iter(i) = mass_iter_i(i);

A(1:N,i) = A_i(1:N);
B(1:N,i) = B_i(1:N);

end

%s

clear Vm_i mass_frac_i mass_iter_i PP_i TT_i
Vm = Old_Vm;


%% Convergence-Velocity
 [Vm,V_flag,V_iter] = converge1(N,M,Vm,Vm_NEW,Vel_tol,V_iter);
Vz = Vm.*cos(phi);

pcolor(z,r,p); shading interp;
colorbar;
hold on
plot(z',r','k');

4


end %Ends while loop for velocity convergance (V_flag)


% %% Radial
% 
% for i=30:M
% 
% [r_NEW_i,VM_NEW_i,RFr_i,Ave_mach_i] = radial(i,N,M,mass_frac,r,Vm,gamma,R,T,delta_z);
% 
% r_NEW(1:N,i) = r_NEW_i(1:N,i);
% Vm_NEW(1:N,i) = VM_NEW_i(1:N,i);
% RFr(i) = RFr_i;
% Ave_mach(i) = Ave_mach_i;
% 
% end
% 
% r_NEW(1:N,1) = r(1:N,1);
% clear r_NEW_i VM_NEW_i RFr_i
% 
% 
% pcolor(z,r_NEW,Vm_NEW); 
% colorbar;
% colormap('jet');
% title('Vm_NEW');
% 
% %% Converge-Radial Change
% [R_flag,Vm,r,r_iter] = converge2(r,r_NEW,M,N,Vm_NEW,Vm,r_iter);
% Vz = Vm.* cos(phi);
% V_tot = sqrt( Vm.^ 2 + V_t.^ 2 );
% 
% end %Ends while loop for radial changes (R_flag)

%% Data Conversion

Vz = Vm.* cos(phi);
Vr = Vm.* sin(phi);
V_tot = sqrt( Vm.^ 2 + V_t.^ 2 );


Pend=0;
Pin=0;
for j = 1:N
    Pend = Pend + p(j,end);
    Pin = Pin + p(j,1);
end

(Pin/N)/100000

(Pend/N)/100000

(Pend/N)/(Pin/N)


pcolor(z,r,p); shading interp;
colorbar;
colormap('jet');
title('p');
hold on
plot(z',r','k');


% xlswrite('Results.xlsx',r,'Radial position')
% xlswrite('Results.xlsx',z,'Axial position')
% xlswrite('Results.xlsx',V_tot,'Total velocity')
% xlswrite('Results.xlsx',rho,'Density')
% xlswrite('Results.xlsx',p,'Pressure')
% xlswrite('Results.xlsx',T,'Temperature')
toc