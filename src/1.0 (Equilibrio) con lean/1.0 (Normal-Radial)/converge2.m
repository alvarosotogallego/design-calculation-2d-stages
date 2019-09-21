function [R_flag,Vm,r,r_iter] =...
converge2(r,r_new,M,N,VM_NEW,Vm,r_iter)
R_flag = 0; %Assume radial changes have converged
r_iter = r_iter + 1;

for i = 2:M
for j = 1:N
X(j,i) = abs(r(j,i) - r_new(j,i));
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
for j = 1:N
Vm(j,i) = VM_NEW(j,i);
r(j,i) = r_new(j,i);
end
end
else
Vm = Vm;
r = r;
end