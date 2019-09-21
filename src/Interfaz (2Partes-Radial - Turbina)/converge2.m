function [R_flag,Vm,r,z,r_iter] =...
converge2(SM,SM_new,M,N,VM_NEW,Vm,r_iter,z,r)
R_flag = 0; %Assume radial changes have converged
r_iter = r_iter + 1;

for i = 2:M
for j = 1:N
X(j,i) = abs(SM(j,i) - SM_new(j,i));
end
end

Xmax = 0; 

for i = 2:M
for j = 1:N
if X(j,i) > Xmax 
    Xmax = X(j,i);
end
end
end

for i = 2:M
for j = 1:N
if X(j,i) >= Xmax & r_iter < 20
R_flag = 1; %Radial Changes have not converged

end
end
end


if R_flag == 1
    
for i = 2:M
    
     Ang = asin((r(N,i)-r(1,i))/SM(N,i));
     Vm(1,i) = VM_NEW(1,i);
     
for j = 2:N
    
Vm(j,i) = VM_NEW(j,i);
SM(j,i) = SM_new(j,i);

r(j,i) = r(j-1,i) + (SM(j,i)-SM(j-1,i))*sin(Ang);
z(j,i) = z(j-1,i) + (SM(j,i)-SM(j-1,i))*cos(Ang);
end
end

end