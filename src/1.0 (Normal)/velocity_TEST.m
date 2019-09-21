function [Vm,mass_frac,mass_iter,Mass_Flow,A,B,V_t,T_0,T,p,p_0,s,H,H_0,V_tot,V_tREL1,V_ang,V_tREL,rho] =...
velocity_TEST(i,SL_Mid,Vm,phi,C,N,dvmdm,dHodq,dsdq,...
drvtdr,r,RFv,rho,mdot_tot,mass_flow_tol,mass_iter,...
FM,FN,T_0,V_t,bladeTest,rel_f,omegaMatriz,Cp,gamma,T,p,p_0,s,H,H_0,R,V_tot,V_tREL1,V_tREL,V_ang,TIPO)

flag_mass = 1;
Vm_MID = Vm(SL_Mid,i);
for j = 1:N
Old_Vm2(j) = Vm(j,i) ^ 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
while flag_mass == 1
VM2NEW(SL_Mid) = Vm_MID ^ 2;

for j=1:N
    if bladeTest(j,i)~=0
        V_t(j,i)=tan(rel_f(j,i))*Vm(j,i)+omegaMatriz(j,i)*r(j,i);
    else
        V_t(j,i)=r(j,i-1)*V_t(j,i-1)/r(j,i);
    end
end

for j = 1:N
A(j) = 2 * (cos(phi(j,i) + psi(j,i)) * C(j,i) + sin(phi(j,i) + ...
psi(j,i)) / Vm(j,i) * dvmdm(j,i));

if TIPO == 2
B(j) = 0.3*2 * (dHodq(j,i) - T_0(j,i) * dsdq(j,i) - 1 / ...
(2 * r(j,i) ^ 2) * drvtdr(j,i) - sin(phi(j,i) + psi(j,i)) * ...
FM(j,i) - cos(phi(j,i) + psi(j,i)) * FN(j,i));
else
    
B(j) = 0.5*2 * (dHodq(j,i) - T_0(j,i) * dsdq(j,i) - 1 / ...
(2 * r(j,i) ^ 2) * drvtdr(j,i) - sin(phi(j,i) + psi(j,i)) * ...
FM(j,i) - cos(phi(j,i) + psi(j,i)) * FN(j,i));
end
end

for j = SL_Mid+1:N
A_Bar(j) = (A(j-1) + A(j)) / 2;
VM2NEW(j) = exp( -A_Bar(j) * (r(j,i) - r(j-1,i)) ) * (VM2NEW(j-1) + ...
((r(j,i) - r(j-1,i)) * B(j-1)) / 2 ) + B(j) * (r(j,i) - ...
r(j-1,i)) / 2;
end
for j = SL_Mid-1:-1:1
A_Bar(j) = (A(j+1) + A(j)) / 2;
VM2NEW(j) = exp( A_Bar(j) * (r(j+1,i) - r(j,i)) ) * (VM2NEW(j+1) + ...
((r(j+1,i) - r(j,i)) * -B(j+1)) / 2) - B(j) * ...
(r(j+1,i) - r(j,i)) / 2;
end

for j = 1:N
VM2NEW(j) = Old_Vm2(j) + RFv * (VM2NEW(j)-Old_Vm2(j));
Old_Vm2(j) = VM2NEW(j);
end

for j = 1:N
VM2NEW(j) = Old_Vm2(j) + RFv * (VM2NEW(j)-Old_Vm2(j));
Old_Vm2(j) = VM2NEW(j);
end
for j = 1:N
%Evitar reflujos
if VM2NEW(j) < 0
VM2NEW(j) = 0.1;
end
%Asignación de la nueva Vm
Vm(j,i) = sqrt(VM2NEW(j));
end
Mass_Flow = 0;
Tube_Flow(1) = 0;

for j=1:N
    if bladeTest(j,i)~=0
        V_tREL(j,i)=tan(rel_f(j,i))*Vm(j,i);
    else
        V_tREL(j,i)=r(j,i-1)*V_t(j,i-1)/r(j,i);
    end
end


for j=1:N
    
    if bladeTest(j,i)~=0
        V_tREL1(j,i)=tan(rel_f(j,i))*Vm(j,i);
    else
        V_tREL1(j,i)=r(j,i-1)*V_t(j,i-1)/r(j,i);
    end
    
    if i < 100
        
        if TIPO == 1
            
    if bladeTest(j,i+1) == 2 & bladeTest(j,i) == 0
        V_tREL1(j,i) = V_tREL1(j,i) - omegaMatriz(j,i+1)*r(j,i+1);
    elseif bladeTest(j,i) == 2 & bladeTest(j,i+1) == 0
        V_tREL1(j,i) = V_tREL1(j,i) + omegaMatriz(j,i)*r(j,i);
    end
    
        else
           
        if bladeTest(j,i+1) == 1 & bladeTest(j,i) == 0
        V_tREL1(j,i) = V_tREL1(j,i) - omegaMatriz(j,i+1)*r(j,i+1);
                
        elseif bladeTest(j,i) == 1 & bladeTest(j,i+1) == 0
        V_tREL1(j,i) = V_tREL1(j,i) + omegaMatriz(j,i)*r(j,i);
        end
    
        end
    
    end
end

for j=1:N
        V_ang(j,i)=omegaMatriz(j,i)*r(j,i);
end

for j=1:N    
    
    if bladeTest(j,i) ~= bladeTest(j,i-1)
        
        T(j,i)= T(j,i-1) - 1/(2*Cp) * ((Vm(j,i)^2 - Vm(j,i-1)^2) + (V_tREL(j,i)^2 - V_tREL1(j,i-1)^2) - ((omegaMatriz(j,i)*r(j,i))^2 - (omegaMatriz(j,i)*r(j,i-1))^2));
        
    else
        
        T(j,i)= T(j,i-1) - 1/(2*Cp) * ((Vm(j,i)^2 - Vm(j,i-1)^2) + (V_tREL(j,i)^2 - V_tREL(j,i-1)^2) - ((omegaMatriz(j,i)*r(j,i))^2 - (omegaMatriz(j,i)*r(j,i-1))^2));
         
    end
        
    V_tot(j,i) = sqrt( Vm(j,i)^2 + V_t(j,i)^2 );

    T_0(j,i) = T(j,i) + (V_tot(j,i)^2) / (2 * Cp);

    p(j,i) = p(j,i-1)*(T(j,i-1)/T(j,i))^(gamma/(1-gamma));

    H_0(j,i) = T_0(j,i)*Cp;

    H(j,i) = H_0(j,i) - (V_tot(j,i)^2)/2; 

    p_0(j,i) = p_0(j,i-1)*(T_0(j,i-1)/T_0(j,i))^(gamma/(1-gamma));
    
    s(j,i) = s(j,i);
    
    rho(j,i) = p(j,i)/ ( R.*T(j,i) );

end


for j = 2:N
Tube_Flow(j) = (rho(j-1,i) + rho(j,i)) / 2 * pi * ((r(j,i)) ^ 2 - ...
(r(j-1,i)) ^ 2) * (Vm(j,i) * cos(phi(j,i)) + Vm(j-1,i) * ...
cos(phi(j-1,i))) / 2;
Mass_Flow = Mass_Flow + Tube_Flow(j);
end

Mass_Flow ;


if abs(Mass_Flow - mdot_tot) < mass_flow_tol
flag_mass = 0; %Flow has converged
mass_frac(1,i) = 0;
mass_frac(N,i) = 1;
for j = 2:N-1
mass_frac(j,i) = Tube_Flow(j) / Mass_Flow + mass_frac(j-1,i);
end
end

2

if Mass_Flow / mdot_tot < 1
    
    4
    
Vm_MID = min([ 1.2 * Vm_MID (2 - Mass_Flow / mdot_tot) * Vm_MID ]);
elseif Mass_Flow / mdot_tot > 1
Vm_MID = max([ 0.8 * Vm_MID (2 - Mass_Flow / mdot_tot) * Vm_MID ]);
end

mass_iter(i) = mass_iter(i) + 1;


end
end %End the while loop for flag_mass