function [Vm,V_flag,V_iter] = converge1(N,M,Vm,Vm_NEW,Vel_tol,V_iter)
V_iter = V_iter + 1;
V_flag = 0; %Assume Velocity has converged, change if not
for i = 2:M
for j = 1:N
X(j,i) = abs(Vm(j,i) - Vm_NEW(j,i));
end
end
for i = 2:M
for j = 1:N
if X(j,i) > Vel_tol
V_flag = 1;
end
end
end
if V_flag == 1
Vm = Vm_NEW;
end
