function [rV_t2,s,V_t,T_0,H,T,p,p_0,V_tot,rho,dsdq,drvtdr,dHodq,dvmdm,...
delta_s,delta_po] = ...
gradient(i,r,T_0,N,Vm,rho,V_t,s,H_0,Cp,R,gamma,T,SL_Mid,...
p_0,V_tot,p,delta_m,rV_t2,M,omega_m,ct,ch,H,test)
for j = SL_Mid:N
omega_r(j,i) = 0;
end
for j = 1:SL_Mid-1;
omega_r(j,i) = 0;
end
for j = 1:N
delta_po1(j,i) = omega_r(j,i)* .5 * V_tot(j,i) ^ 2;
end
%Add in other losses here
delta_po = delta_po1 + 0;
for j = 1:N
if i == 1
delta_s(j,i) = Cp * ((gamma / (gamma-1)) * log((p_0(j,i) - abs(delta_po(j,i))) / p_0(j,i)) - log(T_0(j,i) / T_0(j,i)));
else
delta_s(j,i) = Cp * ((gamma / (gamma-1)) * log((p_0(j,i-1) - abs(delta_po(j,i))) / p_0(j,i-1)) - log(T_0(j,i) / T_0(j,i-1)));
end
end

for j=1:N
if i == 1
s(j,i) = s(j,i) - delta_s(j,i);
else
s(j,i) = s(j,i-1) - delta_s(j,i);
end
end

%Update all other flow field variable according to these new variables
for j = 1:N

V_tot(j,i) = sqrt(Vm(j,i) ^ 2 + V_t(j,i) ^ 2);
H(j,i) = H_0(j,i) - (V_tot(j,i) ^ 2) / 2;
T(j,i) = H(j,i) / Cp;
p(j,i) = T(j,i) ^ (Cp/R) * exp(-s(j,i) / (R));
p_0(j,i) = p(j,i) * (T_0(j,i) / T(j,i)) ^ (gamma / (gamma-1));
rho(j,i) = (p(j,i)) / (R * T(j,i));

end

%Needed in later parts of the program
for j = 1:N
    
    if i < 22
        dHodq(j,i) = 0;
    else
dHodq(j,i) = deriv(j,H_0(1:N,i),r(1:N,i),N,1);
    end
    
end
for j = 1:N
    
    if test == 0
        dsdq(j,i) = 0;
    elseif i < 22
        dsdq(j,i) = 0;
    else
        dsdq(j,i) = deriv(j,s(1:N,i),r(1:N,i),N,1);
    end

end
for j = 1:N
    
        if i < 22
        drvtdr(j,i) = 0;
        else
drvtdr(j,i) = deriv(j,rV_t2(1:N,i),r(1:N,i),N,1);
        end
end
for j = 1:N
if i == 1
dvmdm(j,1) = (delta_m(j,1) + delta_m(j,2)) / delta_m(j,1) * ...
(Vm(j,2)-Vm(j,1))/delta_m(j,2) -(delta_m(j,1) / ...
(delta_m(j,2) + delta_m(j,1))) * (Vm(j,3) - ...
Vm(j,1)) / delta_m(j,2);
elseif i == M
dvmdm(j,M) = (delta_m(j,M) + delta_m(j,M-1)) / delta_m(j,M) * ...
(Vm(j,M)-Vm(j,M-1)) / delta_m(j,M)-(delta_m(j,M) / ...
(delta_m(j,M) + delta_m(j,M-1))) * (Vm(j,M) - ...
Vm(j,M-2)) / delta_m(j,M-1);
else
dvmdm(j,i) = (delta_m(j,i+1) / (delta_m(j,i) + ...
delta_m(j,i+1))) * (Vm(j,i+1) - Vm(j,i)) / ...
delta_m(j,i+1) + (delta_m(j,i+1) / (delta_m(j,i) ...
+delta_m(j,i+1))) * (Vm(j,i) - Vm(j,i-1)) / ...
delta_m(j,i);
end
end
