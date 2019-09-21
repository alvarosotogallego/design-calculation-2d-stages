function plotCircle3D(center,normal,radius,C)

theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
if C == 1
    plot3(points(1,:),points(2,:),points(3,:),'r-');
else
    plot3(points(1,:),points(2,:),points(3,:),'b-');
end
