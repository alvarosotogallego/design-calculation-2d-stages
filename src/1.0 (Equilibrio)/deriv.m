function a = deriv(j,y,x,X,d)
%j is the index of the radial streamline location, used to determine if the derivative is being calculate at an end point or in the middle
%y is the variable for the derivative
%x is locations along the x directions
%X is number of points in the x direction
%d is number of derivative (1st or 2nd)
if d == 1
if j == 1
a = ((x(3) - x(1)) / (x(3) - x(2)) * (y(2) - y(1)) / (x(2) - x(1))) ...
-((x(2) - x(1)) / (x(3) - x(2)) * (y(3) - y(1)) / (x(3) - x(1)));
end
if 1 < j && j < X
a = ((x(j+1)-x(j))/(x(j+1)-x(j-1)))*((y(j)-y(j-1))/(x(j)-x(j-1))) ...
+((x(j)-x(j-1))/(x(j+1)-x(j-1)))*((y(j+1)-y(j))/(x(j+1)-x(j)));
end
if j == X
a = ((x(X)-x(X-2))/(x(X-1)-x(X-2)))*((y(X)-y(X-1))/(x(X)-x(X-1))) ...
-((x(X)-x(X-1))/(x(X-1)-x(X-2)))*((y(X)-y(X-2))/(x(X)-x(X-2)));
end
end
if d == 2
if j == 1
a = 2/(x(3)-x(2))*((y(3)-y(1))/(x(3)-x(1))-(y(2)-y(1))/(x(2)-x(1)));
end
if 1 < j && j < X
a = 2/(x(j+1)-x(j-1))*(((y(j+1)-y(j))/(x(j+1)-x(j)))-...
((y(j)-y(j-1))/(x(j)-x(j-1))));
end
if j == X
a = 2/(x(X-1)-x(X-2))*((y(X)-y(X-1))/(x(X)-x(X-1))-...
(y(X)-y(X-2))/(x(X)-x(X-2)));
end
end
