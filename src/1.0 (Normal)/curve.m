function [C_i,phi_i] = curve(i,r,z,N,M)
if i == M
for j = 1:N
drdz(j,M) = deriv(i,r(j,1:M),z,M,1);
dr2dz2(j,M) = deriv(i,r(j,1:M),z,M,2);
C_i(j,M) = - dr2dz2(j,M) / (( 1 + drdz(j,M) ^ 2) ^ 1.5);
end
elseif i ==1
for j = 1:N
drdz(j,1) = deriv(i,r(j,1:M),z,M,1);
dr2dz2(j,1) = deriv(i,r(j,1:M),z,M,2);
C_i(j,1) = - dr2dz2(j,1) / ((1 + drdz(j,1) ^ 2) ^ 1.5);
end
else
for j = 1:N
dr2dz2(j,i) = deriv(i,r(j,1:M),z,M,2);
drdz(j,i) = deriv(i,r(j,1:M),z,M,1);
C_i(j,i) = - dr2dz2(j,i) / ((1 + drdz(j,i) ^ 2) ^ 1.5);
end
end
for j = 1:N
phi_i(j,i) = atan( drdz(j,i) );
end
end