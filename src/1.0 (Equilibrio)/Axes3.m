function varargout = Axes3(varargin)
% AXES3 MATLAB code for Axes3.fig
%      AXES3, by itself, creates a new AXES3 or raises the existing
%      singleton*.
%
%      H = AXES3 returns the handle to a new AXES3 or the handle to
%      the existing singleton*.
%
%      AXES3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AXES3.M with the given input arguments.
%
%      AXES3('Property','Value',...) creates a new AXES3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Axes3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Axes3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Axes3

% Last Modified by GUIDE v2.5 10-Jun-2019 23:37:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Axes3_OpeningFcn, ...
                   'gui_OutputFcn',  @Axes3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Axes3 is made visible.
function Axes3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Axes3 (see VARARGIN)

gatherAndUpdate(handles)

% Choose default command line output for Axes3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Axes3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Axes3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function gatherAndUpdate(handles)

updateAxes(handles.axes3D);

function updateAxes(axesToUse)
axes(axesToUse)

global PP0 PP1 PP2 PP3 PP4 PP5 PP6 PP7 
global PL1B PL1A PL1BP PL1AP
global PL1B2 PL1A2 PL1B2P PL1A2P 
global PBA1 PBS1 PBA2 PBS2
global P2BA1 P2BS1 P2BA2 P2BS2
global Nalabes

for Ang = 1:90
x = [PP0(1) , PP1(1) , PP2(1) , PP3(1)];
y = [PP0(2) , PP1(2) , PP2(2) , PP3(2)];
z = [ 0 , 0 , 0 , 0];
v = [y;z];
% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta = Ang*4*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

[p1]= [v(1,1) v(2,1);v(1,2) v(2,2);v(1,3) v(2,3);v(1,4) v(2,4)];

n=4;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.02:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL1=l*p1;
PLv3 = [v(3,1);v(3,2);v(3,3);v(3,4)];
PL10=l*PLv3;
plot3(PL1(:,1),PL1(:,2),PL10,'r');
hold on

x = [PP4(1) , PP5(1) , PP6(1) , PP7(1)];
y = [PP4(2) , PP5(2) , PP6(2) , PP7(2)];
z = [ 0 , 0 , 0 , 0];
v = [y;z];
% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta = Ang*4*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

[p1]= [v(1,1) v(2,1);v(1,2) v(2,2);v(1,3) v(2,3);v(1,4) v(2,4)];

n=4;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.02:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL1=l*p1;
PLv3 = [v(3,1);v(3,2);v(3,3);v(3,4)];
PL10=l*PLv3;
plot3(PL1(:,1),PL1(:,2),PL10,'b');
hold on
end

plotCircle3D([PP0(1),0,0],[1,0,0],PP0(2),1)
hold on
plotCircle3D([PP4(1),0,0],[1,0,0],PP4(2),2)
hold on
plotCircle3D([PP3(1),0,0],[1,0,0],PP3(2),1)
hold on
plotCircle3D([PP7(1),0,0],[1,0,0],PP7(2),2)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÁLABE1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAÍZ
for Ang = 1:Nalabes
    
   endPL = length(PL1A);
   endPLB = length(PL1B);
   x = linspace(PBA1(1),PBS1(1),100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A(round(endPL/100*i),2);
   end
   
   for i=1:100
       z(i)=PP0(2);
   end

   central(1)=0;

   for i=2:100
       central(i)=central(i-1)-4*deltax*sin(degtorand*PL1B(round(endPLB/99*(i-1)),2));
   end
   
   for i=1:100
       xdebajo(i)=x(i)-sin(degtorand*PL1B(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centraldebajo(i)=central(i)-0.1*y(i);
   end
   
   for i=1:100
       xarriba(i)=x(i)+sin(degtorand*PL1B(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centralarriba(i)=central(i)+0.1*y(i);
   end
   
    v = [centraldebajo;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
  
   p = plot3(xdebajo,y,z,'Color','r');
   p.LineWidth = 2;
   
   
   hold on
   
   for i=1:100
       z(i)=PP0(2);
   end
   
   v = [centralarriba;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
     

x3_1 = v(1,1);
y3_1 = v(2,1);
z3_1 = v(3,1);

x4_1 = v(1,end);
y4_1 = v(2,end);
z4_1 = v(3,end);

   p = plot3(xarriba,y,z,'Color','r');
   p.LineWidth = 2;
   
  
   hold on

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUNTA
   
   endPL = length(PL1AP);
   endPLB = length(PL1BP);
   x = linspace(PBA2(1),PBS2(1),100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1AP(round(endPL/100*i),2);
   end
   
   for i=1:100
       z(i)=PP4(2);
   end

   central(1)=0;

   for i=2:100
       central(i)=central(i-1)-4*deltax*sin(degtorand*PL1BP(round(endPLB/99*(i-1)),2));
   end
   
   for i=1:100
       xdebajo(i)=x(i)-sin(degtorand*PL1BP(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centraldebajo(i)=central(i)-0.1*y(i);
   end
   
   for i=1:100
       xarriba(i)=x(i)+sin(degtorand*PL1BP(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centralarriba(i)=central(i)+0.1*y(i);
   end
   
    v = [centraldebajo;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
  
   p = plot3(xdebajo,y,z,'Color','b');
   p.LineWidth = 2;
   
   
   hold on
   
   for i=1:100
       z(i)=PP4(2);
   end
   
   v = [centralarriba;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);

x3_2 = v(1,1);
y3_2 = v(2,1);
z3_2 = v(3,1);

x4_2 = v(1,end);
y4_2 = v(2,end);
z4_2 = v(3,end);
     
   p = plot3(xarriba,y,z,'Color','b');
   p.LineWidth = 2;
   hold on
   
       
    px=[x3_1 x3_2];
    py=[y3_1 y3_2];
    pz=[z3_1 z3_2]; 
    
    p = plot3(px,py,pz,'Color','k');
    p.LineWidth = 2;
    hold on
    
    px=[x4_1 x4_2];
    py=[y4_1 y4_2];
    pz=[z4_1 z4_2]; 
    
    p = plot3(px,py,pz,'Color','k');
    p.LineWidth = 2;
    hold on
   
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÁLABE2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAÍZ
for Ang = 1:Nalabes
    
   endPL = length(PL1A2);
   endPLB = length(PL1B2);
   x = linspace(P2BA1(1),P2BS1(1),100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A2(round(endPL/100*i),2);
   end
   
   for i=1:100
       z(i)=PP0(2);
   end

   central(1)=0;

   for i=2:100
       central(i)=central(i-1)-4*deltax*sin(degtorand*PL1B2(round(endPLB/99*(i-1)),2));
   end
   
   for i=1:100
       xdebajo(i)=x(i)-sin(degtorand*PL1B2(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centraldebajo(i)=central(i)-0.1*y(i);
   end
   
   for i=1:100
       xarriba(i)=x(i)+sin(degtorand*PL1B2(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centralarriba(i)=central(i)+0.1*y(i);
   end
   
    v = [centraldebajo;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
  
   p = plot3(xdebajo,y,z,'Color','r');
   p.LineWidth = 2;
   
   
   hold on
   
   for i=1:100
       z(i)=PP0(2);
   end
   
   v = [centralarriba;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
     
x_1 = v(1,1);
y_1 = v(2,1);
z_1 = v(3,1);

x2_1 = v(1,end);
y2_1 = v(2,end);
z2_1 = v(3,end);

   p = plot3(xarriba,y,z,'Color','r');
   p.LineWidth = 2;
   
   
   hold on

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUNTA
   
   endPL = length(PL1A2P);
   endPLB = length(PL1B2P);
   x = linspace(P2BA2(1),P2BS2(1),100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A2P(round(endPL/100*i),2);
   end
   
   for i=1:100
       z(i)=PP4(2);
   end

   central(1)=0;

   for i=2:100
       central(i)=central(i-1)-4*deltax*sin(degtorand*PL1B2P(round(endPLB/99*(i-1)),2));
   end
   
   for i=1:100
       xdebajo(i)=x(i)-sin(degtorand*PL1B2P(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centraldebajo(i)=central(i)-0.1*y(i);
   end
   
   for i=1:100
       xarriba(i)=x(i)+sin(degtorand*PL1B2P(round(endPLB/100*i),2))*0.01;
   end
   for i=1:100
       centralarriba(i)=central(i)+0.1*y(i);
   end
   
    v = [centraldebajo;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
  
   p = plot3(xdebajo,y,z,'Color','b');
   p.LineWidth = 2;
   
   
   hold on
   
   for i=1:100
       z(i)=PP4(2);
   end
   
   v = [centralarriba;z];

% choose a point which will be the center of rotation
y_center = 0;
z_center = 0;
% create a matrix which will be used later in calculations
center = repmat([y_center; z_center], 1, length(x));
% define a x degree counter-clockwise rotation matrix
theta =360/Nalabes*Ang*(pi/180);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
y_rotated = vo(1,:);
z_rotated = vo(2,:);
v = [x;vo];

x = v(1,:);
y = v(2,:);
z = v(3,:);
  
x_2 = v(1,1);
y_2 = v(2,1);
z_2 = v(3,1);

x2_2 = v(1,end);
y2_2 = v(2,end);
z2_2 = v(3,end);

   p = plot3(xarriba,y,z,'Color','b');
   p.LineWidth = 2;
   
   
   hold on

    px=[x_1 x_2];
    py=[y_1 y_2];
    pz=[z_1 z_2]; 
    
    p = plot3(px,py,pz,'Color','k');
    p.LineWidth = 2;
    hold on
   
    px=[x2_1 x2_2];
    py=[y2_1 y2_2];
    pz=[z2_1 z2_2]; 
    
    p = plot3(px,py,pz,'Color','k');
    p.LineWidth = 2;
    hold on
    
end

%%%%%%%%%%%%%%%%%%%%%%% LINEAS ÁLABES





hold off
xlim([-10 (-10+2*PP4(2)+4)])
ylim([(-PP4(2)-3) (PP4(2)+3)])
zlim([(-PP4(2)-3) (PP4(2)+3)])

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
