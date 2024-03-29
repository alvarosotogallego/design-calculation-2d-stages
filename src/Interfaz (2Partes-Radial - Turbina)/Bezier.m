function varargout = Bezier(varargin)
% BEZIER MATLAB code for Bezier.fig
%      BEZIER, by itself, creates a new BEZIER or raises the existing
%      singleton*.
%
%      H = BEZIER returns the handle to a new BEZIER or the handle to
%      the existing singleton*.
%
%      BEZIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEZIER.M with the given input arguments.
%
%      BEZIER('Property','Value',...) creates a new BEZIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bezier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bezier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bezier

% Last Modified by GUIDE v2.5 20-Jun-2019 19:25:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bezier_OpeningFcn, ...
                   'gui_OutputFcn',  @Bezier_OutputFcn, ...
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


% --- Executes just before Bezier is made visible.
function Bezier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bezier (see VARARGIN)

clear global
    
global PP0 PP1 PP2 PP3 PP4 PP5 PP6 PP7 PP8 PP9
global P0 P1 P2 P3 P4 P5 P6 P7 P8 P9
global PLOT1 PLOT2 PLOT3 PLOT4 INLET OUTLET
global CHECK CHECKB TIPO CHECK2 CHECKB2 
global ESPESOR BETA leanangle sweepangle leanangle2 ESPESOR2 BETA2
global gastoM presion0 temperatura0 omega
global PerdidaCT PerdidaCH PerdidaOMEGA
global TIPOAlabe2 TIPOAlabe1 MENUANT2

CHECK = 0;
CHECKB = 0;
ESPESOR = 0;
BETA = 0;
TIPO = 1;
TIPOAlabe1 = 1;
TIPOAlabe2 = 1;
CHECK2 = 0;
CHECKB2 = 0;
ESPESOR2 = 0;
BETA2 = 0;
MENUANT2 = 0;

leanangle =0;
leanangle2 =0;
sweepangle =0;

gastoM = 2.5;      %%  30
presion0 = 20 *100000;   %%  10 *100000; 
temperatura0 = 900;     %%  300
omega = 15000;

PerdidaCT =0;
PerdidaCH =0;
PerdidaOMEGA =0;

axes(handles.axes1)

P0 = impoint(gca,0,5);
setString(P0,'P0');
P1 = impoint(gca,0,1);
setString(P1,'P1');

P2 = impoint(gca,1,1);
setString(P2,'P2');
P3 = impoint(gca,3.2,1);
setString(P3,'P3');        
P4 = impoint(gca,0.7,5);
setString(P4,'P4');
P5 = impoint(gca,0.7,3.5);
setString(P5,'P5');
P6 = impoint(gca,2,3.5);
setString(P6,'P6');
P7 = impoint(gca,3.2,3.5);
setString(P7,'P7');

P8 = impoint(gca,0,6);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUNTOS 8 Y 9
setString(P8,'P8');
P9 = impoint(gca,0.7,6);
setString(P9,'P9');

gatherAndUpdate(handles)

% Choose default command line output for Bezier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bezier wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function gatherAndUpdate(handles)

global PP0 PP1 PP2 PP3 PP4 PP5 PP6 PP7 

gatheredData = gatherData(handles);
updateAxes(handles.axes1,gatheredData);
updateAxes2(handles.axes2,gatheredData);
updateAxes3(handles.axes3,gatheredData);

 set(handles.p0x,'String',num2str(PP0(1)));
 set(handles.p0y,'String',num2str(PP0(2)));
 
 set(handles.p1x,'String',num2str(PP1(1)));
 set(handles.p1y,'String',num2str(PP1(2)));
 
 set(handles.p2x,'String',num2str(PP2(1)));
 set(handles.p2y,'String',num2str(PP2(2)));
 
 set(handles.p3x,'String',num2str(PP3(1)));
 set(handles.p3y,'String',num2str(PP3(2)));
 
 set(handles.p44x,'String',num2str(PP4(1)));
 set(handles.p44y,'String',num2str(PP4(2)));
 
 set(handles.p5x,'String',num2str(PP5(1)));
 set(handles.p5y,'String',num2str(PP5(2)));
 
 set(handles.p6x,'String',num2str(PP6(1)));
 set(handles.p6y,'String',num2str(PP6(2)));
 
 set(handles.p7x,'String',num2str(PP7(1)));
 set(handles.p7y,'String',num2str(PP7(2)));

function gatheredData = gatherData(handles)
gatheredData.SET = get(handles.pushbuttonSET,'Value');
gatheredData.MENU = get(handles.popupmenu1,'Value');
gatheredData.BS1 = get(handles.sliderBS1,'Value');
gatheredData.BS2 = get(handles.sliderBS2,'Value');
gatheredData.BA1 = get(handles.sliderBA1,'Value');
gatheredData.BA2 = get(handles.sliderBA2,'Value');
gatheredData.B2S1 = get(handles.slider2BS1,'Value');
gatheredData.B2S2 = get(handles.slider2BS2,'Value');
gatheredData.B2A1 = get(handles.slider2BA1,'Value');
gatheredData.B2A2 = get(handles.slider2BA2,'Value');


function updateAxes(axesToUse,gd)
axes(axesToUse)

global PP0 PP1 PP2 PP3 PP4 PP5 PP6 PP7 PP8 PP9
global P0 P1 P2 P3 P4 P5 P6 P7 P8 P9
global PLOT1 PLOT2 PLOT3 PLOT4 INLET OUTLET tip hub
global MENUANT MENUANT2
global PBA1 PBA2 PBS1 PBS2
global PLOTPBA1 PLOTPBA2 PLOTPBS1 PLOTPBS2 PLOTBA PLOTBS
global P2BA1 P2BA2 P2BS1 P2BS2
global PLOT2PBA1 PLOT2PBA2 PLOT2PBS1 PLOT2PBS2 PLOT2BA PLOT2BS

global B1 B2 B3 B4 CHECKB PL1B
global PB1 PB2 PB3 PB4
global PLOTB PLOTB2

global B1P B2P B3P B4P PL1BP
global PB1P PB2P PB3P PB4P
global PLOTBP PLOTB2P

global A1 A2 A3 A4 CHECK PL1A TIPO
global PA1 PA2 PA3 PA4
global PLOTA PLOTA2

global A1P A2P A3P A4P PL1AP
global PA1P PA2P PA3P PA4P
global PLOTAP PLOTA2P

global B12 B22 B32 B42 CHECKB2 PL1B2
global PB12 PB22 PB32 PB42
global PLOT2B PLOT2B2

global B12P B22P B32P B42P PL1B2P
global PB12P PB22P PB32P PB42P
global PLOT2BP PLOT2B2P

global A12 A22 A32 A42 CHECK2 PL1A2 
global PA12 PA22 PA32 PA42
global PLOT2A PLOT2A2

global A12P A22P A32P A42P PL1A2P
global PA12P PA22P PA32P PA42P
global PLOT2AP PLOT2A2P


global ESPESOR BETA ESPESOR2 BETA2
global TIPOAlabe2 TIPOAlabe1

global leanangle leanangle2 sweepangle
global r z bladeTest epsilon psi rel_f Ticks02 TicksBeta2 omegaMatriz omega espesorMatriz SM

if gd.MENU == 1
    
    delete(PLOT1)
    delete(PLOT2)
    delete(PLOT3)
    delete(PLOT4)
    delete(INLET)
    delete(OUTLET)
    delete(tip)
    delete(hub)
    
        if MENUANT == 1 & gd.MENU == 1
            P0 = impoint(gca,PP0(1),PP0(2));
            setString(P0,'P0');
            P1 = impoint(gca,PP1(1),PP1(2));
            setString(P1,'P1');
            P2 = impoint(gca,PP2(1),PP2(2));
            setString(P2,'P2');
            P3 = impoint(gca,PP3(1),PP3(2));
            setString(P3,'P3');

            P4 = impoint(gca,PP4(1),PP4(2));
            setString(P4,'P4');
            P5 = impoint(gca,PP5(1),PP5(2));
            setString(P5,'P5');
            P6 = impoint(gca,PP6(1),PP6(2));
            setString(P6,'P6');
            P7 = impoint(gca,PP7(1),PP7(2));
            setString(P7,'P7');
            P8= impoint(gca,PP8(1),PP8(2));
            setString(P8,'P8');
            P9 = impoint(gca,PP9(1),PP9(2));
            setString(P9,'P9');
            
            delete(PLOTPBA1)
            delete(PLOTPBS1)
            delete(PLOTPBA2)
            delete(PLOTPBS2)
            delete(PLOTBA)
            delete(PLOTBS)
            
            delete(PLOT2PBA1)
            delete(PLOT2PBS1)
            delete(PLOT2PBA2)
            delete(PLOT2PBS2)
            delete(PLOT2BA)
            delete(PLOT2BS)
            
        elseif MENUANT == 2 & gd.MENU == 1
            
            P0 = impoint(gca,PP0(1),PP0(2));
            setString(P0,'P0');
            P1 = impoint(gca,PP1(1),PP1(2));
            setString(P1,'P1');
            P2 = impoint(gca,PP2(1),PP2(2));
            setString(P2,'P2');
            P3 = impoint(gca,PP3(1),PP3(2));
            setString(P3,'P3');

            P4 = impoint(gca,PP4(1),PP4(2));
            setString(P4,'P4');
            P5 = impoint(gca,PP5(1),PP5(2));
            setString(P5,'P5');
            P6 = impoint(gca,PP6(1),PP6(2));
            setString(P6,'P6');
            P7 = impoint(gca,PP7(1),PP7(2));
            setString(P7,'P7');
            P8 = impoint(gca,PP8(1),PP8(2));
            setString(P8,'P8');
            P9 = impoint(gca,PP9(1),PP9(2));
            setString(P9,'P9');
            
            delete(PLOTB)
            delete(PLOTB2)
            delete(B1)
            delete(B2)
            delete(B3)
            delete(B4)
            
            delete(PLOTBP)
            delete(PLOTB2P)
            delete(B1P)
            delete(B2P)
            delete(B3P)
            delete(B4P)
            
            CHECKB = 0;
    
        elseif MENUANT == 3 & gd.MENU == 1
         
            P0 = impoint(gca,PP0(1),PP0(2));
            setString(P0,'P0');
            P1 = impoint(gca,PP1(1),PP1(2));
            setString(P1,'P1');
            P2 = impoint(gca,PP2(1),PP2(2));
            setString(P2,'P2');
            P3 = impoint(gca,PP3(1),PP3(2));
            setString(P3,'P3');

            P4 = impoint(gca,PP4(1),PP4(2));
            setString(P4,'P4');
            P5 = impoint(gca,PP5(1),PP5(2));
            setString(P5,'P5');
            P6 = impoint(gca,PP6(1),PP6(2));
            setString(P6,'P6');
            P7 = impoint(gca,PP7(1),PP7(2));
            setString(P7,'P7');
            P8 = impoint(gca,PP8(1),PP8(2));
            setString(P8,'P8');
            P9 = impoint(gca,PP9(1),PP9(2));
            setString(P9,'P9');
            
            delete(PLOTA)
            delete(PLOTA2)
            delete(A1)
            delete(A2)
            delete(A3)
            delete(A4)
            
            delete(PLOTAP)
            delete(PLOTA2P)
            delete(A1P)
            delete(A2P)
            delete(A3P)
            delete(A4P)
            
            
            CHECK = 0;
            
            elseif MENUANT == 4 & gd.MENU == 1
            
            P0 = impoint(gca,PP0(1),PP0(2));
            setString(P0,'P0');
            P1 = impoint(gca,PP1(1),PP1(2));
            setString(P1,'P1');
            P2 = impoint(gca,PP2(1),PP2(2));
            setString(P2,'P2');
            P3 = impoint(gca,PP3(1),PP3(2));
            setString(P3,'P3');

            P4 = impoint(gca,PP4(1),PP4(2));
            setString(P4,'P4');
            P5 = impoint(gca,PP5(1),PP5(2));
            setString(P5,'P5');
            P6 = impoint(gca,PP6(1),PP6(2));
            setString(P6,'P6');
            P7 = impoint(gca,PP7(1),PP7(2));
            setString(P7,'P7');
            P8 = impoint(gca,PP8(1),PP8(2));
            setString(P8,'P8');
            P9 = impoint(gca,PP9(1),PP9(2));
            setString(P9,'P9');
            
            delete(PLOT2B)
            delete(PLOT2B2)
            delete(B12)
            delete(B22)
            delete(B32)
            delete(B42)
            
            delete(PLOT2BP)
            delete(PLOT2B2P)
            delete(B12P)
            delete(B22P)
            delete(B32P)
            delete(B42P)
            
            CHECKB2 = 0;
            
            elseif MENUANT == 5 & gd.MENU == 1
         
            P0 = impoint(gca,PP0(1),PP0(2));
            setString(P0,'P0');
            P1 = impoint(gca,PP1(1),PP1(2));
            setString(P1,'P1');
            P2 = impoint(gca,PP2(1),PP2(2));
            setString(P2,'P2');
            P3 = impoint(gca,PP3(1),PP3(2));
            setString(P3,'P3');

            P4 = impoint(gca,PP4(1),PP4(2));
            setString(P4,'P4');
            P5 = impoint(gca,PP5(1),PP5(2));
            setString(P5,'P5');
            P6 = impoint(gca,PP6(1),PP6(2));
            setString(P6,'P6');
            P7 = impoint(gca,PP7(1),PP7(2));
            setString(P7,'P7');
            P8 = impoint(gca,PP8(1),PP8(2));
            setString(P8,'P8');
            P9 = impoint(gca,PP9(1),PP9(2));
            setString(P9,'P9');
            
            delete(PLOT2A)
            delete(PLOT2A2)
            delete(A12)
            delete(A22)
            delete(A32)
            delete(A42)
            
            delete(PLOT2AP)
            delete(PLOT2A2P)
            delete(A12P)
            delete(A22P)
            delete(A32P)
            delete(A42P)
            
            CHECK2 = 0;
            
        end
        
    
    PP0 = getPosition(P0);
    PP1 = getPosition(P1);
    PP2 = getPosition(P2);
    PP3 = getPosition(P3);
    PP4 = getPosition(P4);
    PP5 = getPosition(P5);
    PP6 = getPosition(P6);
    PP7 = getPosition(P7);
    PP8 = getPosition(P8);
    PP9 = getPosition(P9);

[p1]= [PP0(1) PP0(2);PP1(1) PP1(2);PP2(1) PP2(2);PP3(1) PP3(2)];
[p2]= [PP4(1) PP4(2);PP5(1) PP5(2);PP6(1) PP6(2);PP7(1) PP7(2)];
[p3]= [PP8(1) PP8(2);PP0(1) PP0(2)];
[p4]= [PP9(1) PP9(2);PP4(1) PP4(2)];

n=4;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.0002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL1=l*p1;
PLOT1 = line(PL1(:,1),PL1(:,2));
hold on
PLOT2 = line(p1(:,1),p1(:,2));
hold on
PL2=l*p2;
PLOT3 = line(PL2(:,1),PL2(:,2));
hold on
PLOT4 = line(p2(:,1),p2(:,2));
hold on
INLET = plot([PP8(1) PP9(1)],[PP8(2) PP9(2)]);
hold on
OUTLET = plot([PP3(1) PP7(1)],[PP3(2) PP7(2)]);
hold on
tip = plot([PP8(1) PP0(1)],[PP8(2) PP0(2)]);
hold on
hub = plot([PP9(1) PP4(1)],[PP9(2) PP4(2)]);

n=2;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL3=l*p3;
PL4=l*p4;

Lgth1 = length(PL1);
Lgth2 = length(PL2);
Lgth3 = length(PL3);
Lgth4 = length(PL4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=zeros([50 100]);

z1=zeros([50 20]);
z2=zeros([50 80]);

for j = 1:80
    
    if j==1
    inter=linspace(PL1(1,1), PL2(1,1), 50);
    elseif j==80
    inter=linspace(PL1(end,1), PL2(end,1), 50);
    else
    inter=linspace(PL1(round(-1/79*Lgth1+1/79*j*Lgth1),1), PL2(round(-1/79*Lgth2+1/79*j*Lgth2),1), 50);
    end
    
       for i = 1:50
           z2(i,j) = inter(i);
       end
       
end

for j = 1:20 
    
    if j==1
    inter=linspace(PL3(1,1), PL4(1,1), 50);
    elseif j==20
    inter=linspace(PL3(end,1), PL4(end,1), 50);
    else
    inter=linspace(PL3(round(-1/19*Lgth3+1/19*j*Lgth3),1), PL4(round(-1/19*Lgth4+1/19*j*Lgth4),1), 50);
    end
    
       for i = 1:50
           z1(i,j) = inter(i);
       end
       
end

z = [z1 z2];
z = z*0.01;

r=zeros([50 100]);

r1=zeros([50 20]);
r2=zeros([50 80]);

for j = 1:80
    
    if j==1
    inter=linspace(PL1(1,2), PL2(1,2), 50);
    elseif j==80
    inter=linspace(PL1(end,2), PL2(end,2), 50);
    else
    inter=linspace(PL1(round(-1/79*Lgth1+1/79*j*Lgth1),2), PL2(round(-1/79*Lgth2+1/79*j*Lgth2),2), 50);
    end
    
       for i = 1:50
           r2(i,j) = inter(i);
       end
       
end

for j = 1:20
    
    if j==1
    inter=linspace(PL3(1,2), PL4(1,2), 50);
    elseif j==20
    inter=linspace(PL3(end,2), PL4(end,2), 50);
    else
    inter=linspace(PL3(round(-1/19*Lgth3+1/19*j*Lgth3),2), PL4(round(-1/19*Lgth4+1/19*j*Lgth4),2), 50);
    end
    
       for i = 1:50
           r1(i,j) = inter(i);
       end
       
end

r = [r1 r2];
r=r*0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis([-2 8 -2 8])

elseif gd.MENU == 2
    
    CHECK = 0;
    CHECKB = 0;
    delete(PLOT2)
    delete(PLOT4)
    delete(P0)
    delete(P1)
    delete(P2)
    delete(P3)
    delete(P4)
    delete(P5)
    delete(P6)
    delete(P7)
    delete(P8)
    delete(P9)
    delete(PLOTPBA1)
    delete(PLOTPBS1)
    delete(PLOTPBA2)
    delete(PLOTPBS2)
    delete(PLOTBA)
    delete(PLOTBS)
    delete(PLOT2PBA1)
    delete(PLOT2PBS1)
    delete(PLOT2PBA2)
    delete(PLOT2PBS2)
    delete(PLOT2BA)
    delete(PLOT2BS)
    
    delete(PLOTA)
    delete(PLOTA2)
    delete(A1)
    delete(A2)
    delete(A3)
    delete(A4)
    
    delete(PLOTAP)
    delete(PLOTA2P)
    delete(A1P)
    delete(A2P)
    delete(A3P)
    delete(A4P)
    
    delete(PLOTB)
    delete(PLOTB2)
    delete(B1)
    delete(B2)
    delete(B3)
    delete(B4)
    
    delete(PLOTBP)
    delete(PLOTB2P)
    delete(B1P)
    delete(B2P)
    delete(B3P)
    delete(B4P)
    
    delete(PLOT2A)
    delete(PLOT2A2)
    delete(A12)
    delete(A22)
    delete(A32)
    delete(A42)
    
    delete(PLOT2AP)
    delete(PLOT2A2P)
    delete(A12P)
    delete(A22P)
    delete(A32P)
    delete(A42P)
    
    delete(PLOT2B)
    delete(PLOT2B2)
    delete(B12)
    delete(B22)
    delete(B32)
    delete(B42)
    
    delete(PLOT2BP)
    delete(PLOT2B2P)
    delete(B12P)
    delete(B22P)
    delete(B32P)
    delete(B42P)
    
    [p1]= [PP0(1) PP0(2);PP1(1) PP1(2);PP2(1) PP2(2);PP3(1) PP3(2)];
    [p2]= [PP4(1) PP4(2);PP5(1) PP5(2);PP6(1) PP6(2);PP7(1) PP7(2)];

    n=4;
    n1=n-1;

    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1=l*p1;
    PL2=l*p2;
    
    Lgth1 = length(PL1);
    Lgth2 = length(PL2);
    
[p3]= [PP8(1) PP8(2);PP0(1) PP0(2)];
[p4]= [PP9(1) PP9(2);PP4(1) PP4(2)];

n=2;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL3=l*p3;
PL4=l*p4;

Lgth3 = length(PL3);
Lgth4 = length(PL4);
    
    PBA1 = [PL3(round(gd.BA1*Lgth3),1) PL3(round(gd.BA1*Lgth3),2)];
    PBS1 = [PL3(round(gd.BS1*Lgth3),1) PL3(round(gd.BS1*Lgth3),2)];
    PBA2 = [PL4(round(gd.BA2*Lgth4),1) PL4(round(gd.BA2*Lgth4),2)];
    PBS2 = [PL4(round(gd.BS2*Lgth4),1) PL4(round(gd.BS2*Lgth4),2)];
    
    P2BA1 = [PL1(round(gd.B2A1*Lgth1),1) PL1(round(gd.B2A1*Lgth1),2)];
    P2BS1 = [PL1(round(gd.B2S1*Lgth1),1) PL1(round(gd.B2S1*Lgth1),2)];
    P2BA2 = [PL2(round(gd.B2A2*Lgth2),1) PL2(round(gd.B2A2*Lgth2),2)];
    P2BS2 = [PL2(round(gd.B2S2*Lgth2),1) PL2(round(gd.B2S2*Lgth2),2)];
    
    hold on
    PLOTPBS1 = plot(PBS1(1),PBS1(2),'r*');
    hold on
    PLOTPBS2 = plot(PBS2(1),PBS2(2),'r*');
    hold on
    PLOTPBA1 = plot(PBA1(1),PBA1(2),'b*');
    hold on
    PLOTPBA2 = plot(PBA2(1),PBA2(2),'b*');
    hold on
    PLOTBA = plot([PBA1(1) PBA2(1)],[PBA1(2) PBA2(2)]);
    hold on
    PLOTBS = plot([PBS1(1) PBS2(1)],[PBS1(2) PBS2(2)]);
    hold on
    PLOT2PBS1 = plot(P2BS1(1),P2BS1(2),'r*');
    hold on
    PLOT2PBS2 = plot(P2BS2(1),P2BS2(2),'r*');
    hold on
    PLOT2PBA1 = plot(P2BA1(1),P2BA1(2),'b*');
    hold on
    PLOT2PBA2 = plot(P2BA2(1),P2BA2(2),'b*');
    hold on
    PLOT2BA = plot([P2BA1(1) P2BA2(1)],[P2BA1(2) P2BA2(2)]);
    hold on
    PLOT2BS = plot([P2BS1(1) P2BS2(1)],[P2BS1(2) P2BS2(2)]);

    if MENUANT2 == 2 || MENUANT2 == 3 || MENUANT2 == 4 || MENUANT2 == 5

[p1]= [PP0(1) PP0(2);PP1(1) PP1(2);PP2(1) PP2(2);PP3(1) PP3(2)];
[p2]= [PP4(1) PP4(2);PP5(1) PP5(2);PP6(1) PP6(2);PP7(1) PP7(2)];

n=4;
n1=n-1;

for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];

for u=0:0.0002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end

PL1=l*p1;
PLOT1 = line(PL1(:,1),PL1(:,2));
hold on
PL2=l*p2;
PLOT3 = line(PL2(:,1),PL2(:,2));
hold on
INLET = plot([PP8(1) PP9(1)],[PP8(2) PP9(2)]);
hold on
OUTLET = plot([PP3(1) PP7(1)],[PP3(2) PP7(2)]);
hold on
tip = plot([PP8(1) PP0(1)],[PP8(2) PP0(2)]);
hold on
hub = plot([PP9(1) PP4(1)],[PP9(2) PP4(2)]);
        
    end
   
% axis([PP0(1)-3 PP7(1)+3 PP0(2)-3 PP7(2)+3])
axis([-2 8 -2 8])
    
    if ESPESOR == 1 & BETA == 1 & BETA2 == 1 & ESPESOR2 == 1
        
    %Creaci�n de la matriz de datos del �labe
    
    bladeTest=zeros([50 100]);
    epsilon=zeros([50 100]);
    psi=zeros([50 100]);
    rel_f=zeros([50 100]);
    omegaMatriz=zeros([50 100]);
    espesorMatriz=zeros([50 100]);
    
        %Borde de ataque1
    
    for i = 1:50
    TicksBeta(i) = 0;
    Ticks0(i) = 0;
    TicksBeta2(i) = 0;
    Ticks02(i) = 0;
    end
    
    if PBA1(1)==PBA2(1)
        for i = 1:50
            xq(i) = PBA1(1);
        end
        4
        for i = 1:50
            for j = 1:100
                if z(i,j)*100 < xq(i) 
                    bladeTest(i,j) = 0;
                    Ticks0(i) = Ticks0(i) + 1;
                else
                    bladeTest(i,j) = 1;
                    epsilon(i,j) = leanangle;
                    psi(i,j) = sweepangle;
                    TicksBeta(i) = TicksBeta(i) + 1;
                    
                    if TIPO == 2
                      omegaMatriz(i,j)=omega*2*3.1416/60;
                    end
                        
                end
            end
        end
        
    else
        
        x = [PBA1(1) PBA2(1)];
        y = [PBA1(2) PBA2(2)];
        xq = linspace(PBA1(2),PBA2(2),50);
        
        for i = 1:50
            for j = 1:100
                if r(i,j)*100 > xq(i) 
                    bladeTest(i,j) = 0;
                    Ticks0(i) = Ticks0(i) + 1;
                else
                    bladeTest(i,j) = 1;
                    epsilon(i,j) = leanangle;
                    psi(i,j) = sweepangle;
                    TicksBeta(i) = TicksBeta(i) + 1;
                    
                    if TIPO == 2
                      omegaMatriz(i,j)=omega*2*3.1416/60;
                    end
                    
                end
            end
        end
    end
    
    %Borde de salida1
   
    if PBS1(1)==PBS2(1)
        for i = 1:50
            xq(i) = PBS1(1);
        end
        4
        for i = 1:50
            for j = 100:-1:1
                if z(i,j)*100 < xq(i)
                    bladeTest(i,j) = 0;
                    epsilon(i,j) = 0;
                    psi(i,j) = 0;
                    TicksBeta(i) = TicksBeta(i) - 1;
                    if TIPO == 2
                    omegaMatriz(i,j)=0;
                    end
                end
            end
        end
        
    else
        x = [PBS1(1) PBS2(1)];
        y = [PBS1(2) PBS2(2)];
        xq = linspace(PBS1(2),PBS2(2),50);

        for i = 1:50
            for j = 100:-1:1
                if r(i,j)*100 < xq(i) 
                    bladeTest(i,j) = 0;
                    epsilon(i,j) = 0;
                    psi(i,j) = 0;
                    TicksBeta(i) = TicksBeta(i) - 1;
                    if TIPO == 2
                    omegaMatriz(i,j)=0;
                    end
                end
            end
        end
        
    end

    
    for i = 1:50 %%%%%%%%%%%%%% ASIGNACI�N BETA �LABE 1
        n = TicksBeta(i);
        endPLB = length(PL1B);
        degtorand = 3.14166/180;
        
        for k = 1:endPLB
        PL1Belse(k)=round(i*(PL1BP(k,2)-PL1B(k,2))/50);
        end
        
        if i == 1
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1B(1,2);
            elseif j == n
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1B(end,2);
            else
               rel_f(i,j+Ticks0(i)) = degtorand*PL1B(round(endPLB/n*j),2);
            end
        end
        elseif i == 50
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1BP(1,2);
            elseif j == n
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1BP(end,2);
            else
               rel_f(i,j+Ticks0(i)) = degtorand*PL1BP(round(endPLB/n*j),2);
            end
        end
        else
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1B(1,2)+degtorand*PL1Belse(1);
            elseif j == n
               rel_f(i,j+Ticks0(i)) =  degtorand*PL1B(end,2)+degtorand*PL1Belse(end);
            else
               rel_f(i,j+Ticks0(i)) = degtorand*PL1B(round(endPLB/n*j),2)+degtorand*PL1Belse(round(endPLB/n*j));
            end
        end   
        end
        
    end
    
    
    for i = 1:50 %%%%%%%%%%%%%% ASIGNACI�N BETA ESPESOR 1
        n = TicksBeta(i);
        endPLA = length(PL1A);
        
        for k = 1:endPLA
        PL1Aelse(k)=round(i*(PL1AP(k,2)-PL1A(k,2))/50,2);
        end
        
        if i == 1
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks0(i)) =  PL1A(1,2);
            elseif j == n
               espesorMatriz(i,j+Ticks0(i)) =  PL1A(end,2);
            else
               espesorMatriz(i,j+Ticks0(i)) = PL1A(round(endPLA/n*j),2);
            end
        end
        elseif i == 50
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks0(i)) =  PL1AP(1,2);
            elseif j == n
               espesorMatriz(i,j+Ticks0(i)) =  PL1AP(end,2);
            else
               espesorMatriz(i,j+Ticks0(i)) = PL1AP(round(endPLA/n*j),2);
            end
        end
        else
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks0(i)) =  PL1A(1,2)+PL1Aelse(1);
            elseif j == n
               espesorMatriz(i,j+Ticks0(i)) =  PL1A(end,2)+PL1Aelse(end);
            else
               espesorMatriz(i,j+Ticks0(i)) = PL1A(round(endPLA/n*j),2)+PL1Aelse(i)*j/100;
            end
        end   
        end
        
    end
    
    
            i_alabe = TicksBeta(1)+Ticks0(1);
    
    %Borde de ataque2
    
    if P2BA1(1)==P2BA2(1)
        for i = 1:50
            xq(i) = P2BA1(1);
        end
        
        for i = 1:50
            for j = 1:100
                if z(i,j)*100 < xq(i) 
                    Ticks02(i) = Ticks02(i) + 1;
                else
                    bladeTest(i,j) = 2;
                    
                    if TIPO == 1
                    omegaMatriz(i,j)=omega*2*3.1416/60;
                    end
                    
                    epsilon(i,j) = leanangle2;
                    psi(i,j) = sweepangle;
                    TicksBeta2(i) = TicksBeta2(i) + 1;
                end
            end
        end
        
    else
        
        x = [P2BA1(1) P2BA2(1)];
        y = [P2BA1(2) P2BA2(2)];
        xq = linspace(P2BA1(1),P2BA2(1),50);
        
        for i = 1:50
            for j = 1:100
                if z(i,j)*100 < xq(i) 
                    Ticks02(i) = Ticks02(i) + 1;
                else
                    bladeTest(i,j) = 2;
                    
                    if TIPO == 1
                    omegaMatriz(i,j)=omega*2*3.1416/60;
                    end
                    
                    epsilon(i,j) = leanangle2;
                    psi(i,j) = sweepangle;
                    TicksBeta2(i) = TicksBeta2(i) + 1;
                end
            end
        end
    end
    
        %Borde de salida2
   
    if P2BS1(1)==P2BS2(1)
        for i = 1:50
            xq(i) = P2BS1(1);
        end
        
        for i = 1:50
            for j = 100:-1:1
                if z(i,j)*100 > xq(i)
                    bladeTest(i,j) = 0;
                    
                    if TIPO == 1
                    omegaMatriz(i,j)=0;
                    end
                    
                    epsilon(i,j) = 0;
                    psi(i,j) = 0;
                    TicksBeta2(i) = TicksBeta2(i) - 1;
                end
            end
        end
        
    else
        x = [P2BS1(1) P2BS2(1)];
        y = [P2BS1(2) P2BS2(2)];
        xq = linspace(P2BS1(1),P2BS2(1),50);

        for i = 1:50
            for j = 100:-1:1
                if z(i,j)*100 > xq(i) 
                    bladeTest(i,j) = 0;
                    
                    if TIPO == 1
                    omegaMatriz(i,j)=0;
                    end
                    
                    epsilon(i,j) = 0;
                    psi(i,j) = 0;
                    TicksBeta2(i) = TicksBeta2(i) - 1;
                end
            end
        end
        
    end

    for i = 1:50 %%%%%%%%%%%%%% ASIGNACI�N BETA �LABE 2
        n = TicksBeta2(i);
        endPLB = length(PL1B2);
        degtorand = 3.14166/180;
        
        for k = 1:endPLB
        PL1Belse(k)=round(i*(PL1B2P(k,2)-PL1B2(k,2))/50);
        end
        
        if i == 1
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2(1,2);
            elseif j == n
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2(end,2);
            else
               rel_f(i,j+Ticks02(i)) = degtorand*PL1B2(round(endPLB/n*j),2);
            end
        end
        elseif i == 50
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2P(1,2);
            elseif j == n
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2P(end,2);
            else
               rel_f(i,j+Ticks02(i)) = degtorand*PL1B2P(round(endPLB/n*j),2);
            end
        end
        else
        for j = 1:n  
            if  j == 1              
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2(1,2)+degtorand*PL1Belse(1);
            elseif j == n
               rel_f(i,j+Ticks02(i)) =  degtorand*PL1B2(end,2)+degtorand*PL1Belse(end);
            else
               rel_f(i,j+Ticks02(i)) = degtorand*PL1B2(round(endPLB/n*j),2)+degtorand*PL1Belse(round(endPLB/n*j));
            end
        end   
        end
        
    end
        
   
     for i = 1:50 %%%%%%%%%%%%%% ASIGNACI�N BETA ESPESOR 2
        n = TicksBeta2(i);
        endPLA = length(PL1A2);
        
        for k = 1:endPLA
        PL1Aelse(k)=round(i*(PL1A2P(k,2)-PL1A2(k,2))/50,2);
        end
        
        if i == 1
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2(1,2);
            elseif j == n
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2(end,2);
            else
               espesorMatriz(i,j+Ticks02(i)) = PL1A2(round(endPLA/n*j),2);
            end
        end
        elseif i == 50
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2P(1,2);
            elseif j == n
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2P(end,2);
            else
               espesorMatriz(i,j+Ticks02(i)) = PL1A2P(round(endPLA/n*j),2);
            end
        end
        else
        for j = 1:n  
            if  j == 1              
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2(1,2)+PL1Aelse(1);
            elseif j == n
               espesorMatriz(i,j+Ticks02(i)) =  PL1A2(end,2)+PL1Aelse(end);
            else
               espesorMatriz(i,j+Ticks02(i)) = PL1A2(round(endPLA/n*j),2)++PL1Aelse(i)*j/100;
            end
        end   
        end
        
    end
    
     espesorMatriz = espesorMatriz/1000;
     
          i_alabe2 = Ticks02(1);
          i_alabe3 = Ticks02(1)+TicksBeta2(1); 
          
  
          
            %%%%%%%%%%%%%%%%% GEOMETR�A m y tetha
          
          
          SM   = zeros([50 100]);
          
          for i = 1:100
          for j = 2:50
              SM(j,i)  = SM(j-1,i) + sqrt((r(j,i)-r(j-1,i))^2+(z(j,i)-z(j-1,i))^2);
          end
          end

          
          %%%%%%%%%%%%%%%%%
          
          
    end
    
    
    
    
elseif gd.MENU == 3 %%%%%%%%%%%% BETA �LABE 1
    
    CHECK = 0;
    delete(PLOT1)
    delete(PLOT2)
    delete(PLOT3)
    delete(PLOT4)
    delete(INLET)
    delete(OUTLET)
    delete(P0)
    delete(P1)
    delete(P2)
    delete(P3)
    delete(P4)
    delete(P5)
    delete(P6)
    delete(P7)
    delete(P8)
    delete(P9)
    delete(PLOTPBA1)
    delete(PLOTPBS1)
    delete(PLOTPBA2)
    delete(PLOTPBS2)
    delete(PLOTBA)
    delete(PLOTBS)
    delete(PLOT2PBA1)
    delete(PLOT2PBS1)
    delete(PLOT2PBA2)
    delete(PLOT2PBS2)
    delete(PLOT2BA)
    delete(PLOT2BS)
    
    delete(PLOTA)
    delete(PLOTA2)
    delete(A1)
    delete(A2)
    delete(A3)
    delete(A4)
    
    delete(PLOTAP)
    delete(PLOTA2P)
    delete(A1P)
    delete(A2P)
    delete(A3P)
    delete(A4P)
    
    delete(PLOTB)
    delete(PLOTB2)
    
    delete(PLOTBP)
    delete(PLOTB2P)
    
    delete(PLOT2A)
    delete(PLOT2A2)
    delete(A12)
    delete(A22)
    delete(A32)
    delete(A42)
    
    delete(PLOT2AP)
    delete(PLOT2A2P)
    delete(A12P)
    delete(A22P)
    delete(A32P)
    delete(A42P)
    
    delete(PLOT2B)
    delete(PLOT2B2)
    delete(B12)
    delete(B22)
    delete(B32)
    delete(B42)
    
    delete(PLOT2BP)
    delete(PLOT2B2P)
    delete(B12P)
    delete(B22P)
    delete(B32P)
    delete(B42P)
    
    if CHECKB == 1
        
        PB1 = getPosition(B1);   %%%%%%%%%%%% RA�Z
        PB2 = getPosition(B2);
        PB3 = getPosition(B3);
        PB4 = getPosition(B4);
       
        delete(B1)
        delete(B2)
        delete(B3)
        delete(B4)
        
        B1 = impoint(gca,0,PB1(2));
        setString(B1,'R1');
        B2 = impoint(gca,PB2(1),PB2(2));
        setString(B2,'R2');
        B3 = impoint(gca,PB3(1),PB3(2));
        setString(B3,'R3');
        B4 = impoint(gca,1,PB4(2));
        setString(B4,'R4');
        
        PB1P = getPosition(B1P);   %%%%%%%%%%%% PUNTA
        PB2P = getPosition(B2P);
        PB3P = getPosition(B3P);
        PB4P = getPosition(B4P);
       
        delete(B1P)
        delete(B2P)
        delete(B3P)
        delete(B4P)
        
        B1P = impoint(gca,0,PB1P(2));
        setString(B1P,'P1');
        B2P = impoint(gca,PB2P(1),PB2P(2));
        setString(B2P,'P2');
        B3P = impoint(gca,PB3P(1),PB3P(2));
        setString(B3P,'P3');
        B4P = impoint(gca,1,PB4P(2));
        setString(B4P,'P4');
        
    elseif CHECKB == 0
        
        
        if TIPO == 1
        
        B1 = impoint(gca,0,0);%%%%%%%% RA�Z
        setString(B1,'R1');
        B2 = impoint(gca,0.25,0);
        setString(B2,'R2');
        B3 = impoint(gca,0.5,0);
        setString(B3,'R3');
        B4 = impoint(gca,1,0);
        setString(B4,'R4');
        
        B1P = impoint(gca,0,0);%%%%%%%% PUNTA
        setString(B1P,'P1');
        B2P = impoint(gca,0.25,0);
        setString(B2P,'P2');
        B3P = impoint(gca,0.5,0);
        setString(B3P,'P3');
        B4P = impoint(gca,1,0);
        setString(B4P,'P4');
        
        else
            
        B1 = impoint(gca,0,-20);%%%%%%%%%%% RA�Z
        setString(B1,'R1');
        B2 = impoint(gca,0.25,-10);
        setString(B2,'R2');
        B3 = impoint(gca,0.5,0);
        setString(B3,'R3');
        B4 = impoint(gca,1,20);
        setString(B4,'R4');
        
        B1P = impoint(gca,0,-20);%%%%%%%%%%% PUNTA
        setString(B1P,'P1');
        B2P = impoint(gca,0.25,-15);
        setString(B2P,'P2');
        B3P = impoint(gca,0.5,-10);
        setString(B3P,'P3');
        B4P = impoint(gca,1,0);
        setString(B4P,'P4');
        
        end
        
        PB1 = getPosition(B1);%%%%%%%%%%% RA�Z
        PB2 = getPosition(B2);
        PB3 = getPosition(B3);
        PB4 = getPosition(B4);
        
        PB1P = getPosition(B1P);%%%%%%%%%%% PUNTA
        PB2P = getPosition(B2P);
        PB3P = getPosition(B3P);
        PB4P = getPosition(B4P);
        
        CHECKB = 1;
        
    end
    
        [p1]= [0 PB1(2);PB2(1) PB2(2);PB3(1) PB3(2);1 PB4(2)];%%%%%%%%%%% RA�Z
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1B=l*p1;
    PLOTB = line(PL1B(:,1),PL1B(:,2));
    hold on
    PLOTB2 = line(p1(:,1),p1(:,2));
    
    
            [p1]= [0 PB1P(2);PB2P(1) PB2P(2);PB3P(1) PB3P(2);1 PB4P(2)];%%%%%%%%%%% PUNTA
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1BP=l*p1;
    hold on
    PLOTBP = line(PL1BP(:,1),PL1BP(:,2));
    hold on
    PLOTB2P = line(p1(:,1),p1(:,2));
    
    axis([-0.25 1.25 -60 90])
    
elseif gd.MENU == 4   %%%%%%% ESPESOR �LABE 1
    
    CHECKB = 0;
    delete(PLOT1)
    delete(PLOT2)
    delete(PLOT3)
    delete(PLOT4)
    delete(INLET)
    delete(OUTLET)
    delete(P0)
    delete(P1)
    delete(P2)
    delete(P3)
    delete(P4)
    delete(P5)
    delete(P6)
    delete(P7)
    delete(P8)
    delete(P9)
    delete(PLOTPBA1)
    delete(PLOTPBS1)
    delete(PLOTPBA2)
    delete(PLOTPBS2)
    delete(PLOTBA)
    delete(PLOTBS)
    delete(PLOT2PBA1)
    delete(PLOT2PBS1)
    delete(PLOT2PBA2)
    delete(PLOT2PBS2)
    delete(PLOT2BA)
    delete(PLOT2BS)
    
    delete(PLOTA)
    delete(PLOTA2)
    
    delete(PLOTAP)
    delete(PLOTA2P)
    
    delete(PLOTB)
    delete(PLOTB2)
    delete(B1)
    delete(B2)
    delete(B3)
    delete(B4)
    
    delete(PLOTBP)
    delete(PLOTB2P)
    delete(B1P)
    delete(B2P)
    delete(B3P)
    delete(B4P)
    
    delete(PLOT2A)
    delete(PLOT2A2)
    delete(A12)
    delete(A22)
    delete(A32)
    delete(A42)
    
    delete(PLOT2AP)
    delete(PLOT2A2P)
    delete(A12P)
    delete(A22P)
    delete(A32P)
    delete(A42P)
    
    delete(PLOT2B)
    delete(PLOT2B2)
    delete(B12)
    delete(B22)
    delete(B32)
    delete(B42)
    
    delete(PLOT2BP)
    delete(PLOT2B2P)
    delete(B12P)
    delete(B22P)
    delete(B32P)
    delete(B42P)
    
    if CHECK == 1
        
        PA1 = getPosition(A1);%%%%%%%%%%%  RA�Z
        PA2 = getPosition(A2);
        PA3 = getPosition(A3);
        PA4 = getPosition(A4);
        
        delete(A1)
        delete(A2)
        delete(A3)
        delete(A4)
        
        A1 = impoint(gca,0,PA1(2));
        setString(A1,'R1');
        A2 = impoint(gca,PA2(1),PA2(2));
        setString(A2,'R2');
        A3 = impoint(gca,PA3(1),PA3(2));
        setString(A3,'R3');
        A4 = impoint(gca,1,PA4(2));
        setString(A4,'R4');
        
        PA1P = getPosition(A1P);%%%%%%%%%%%  PUNTA
        PA2P = getPosition(A2P);
        PA3P = getPosition(A3P);
        PA4P = getPosition(A4P);
        
        delete(A1P)
        delete(A2P)
        delete(A3P)
        delete(A4P)
        
        A1P = impoint(gca,0,PA1P(2));
        setString(A1P,'P1');
        A2P = impoint(gca,PA2P(1),PA2P(2));
        setString(A2P,'P2');
        A3P = impoint(gca,PA3P(1),PA3P(2));
        setString(A3P,'P3');
        A4P = impoint(gca,1,PA4P(2));
        setString(A4P,'P4');
        
    elseif CHECK == 0
        
        A1 = impoint(gca,0,2);%%%%%%%%%%%  RA�Z
        setString(A1,'R1');
        A2 = impoint(gca,0.2,6);
        setString(A2,'R2');
        A3 = impoint(gca,0.6,2);
        setString(A3,'R3');
        A4 = impoint(gca,1,1);
        setString(A4,'R4');
        
        PA1 = getPosition(A1);
        PA2 = getPosition(A2);
        PA3 = getPosition(A3);
        PA4 = getPosition(A4);
        
        A1P = impoint(gca,0,2);%%%%%%%%%%%  PUNTA
        setString(A1P,'P1');
        A2P = impoint(gca,0.2,6);
        setString(A2P,'P2');
        A3P = impoint(gca,0.6,2);
        setString(A3P,'P3');
        A4P = impoint(gca,1,1);
        setString(A4P,'P4');
        
        PA1P = getPosition(A1P);
        PA2P = getPosition(A2P);
        PA3P = getPosition(A3P);
        PA4P = getPosition(A4P);
        
        CHECK = 1;
        
    end
    
    
    [p1]= [0 PA1(2);PA2(1) PA2(2);PA3(1) PA3(2);1 PA4(2)];%%%%%%%%%%%  RA�Z
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1A=l*p1;
    PLOTA = line(PL1A(:,1),PL1A(:,2));
    hold on
    PLOTA2 = line(p1(:,1),p1(:,2));
    
        [p1]= [0 PA1P(2);PA2P(1) PA2P(2);PA3P(1) PA3P(2);1 PA4P(2)];%%%%%%%%%%%  PUNTA
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1AP=l*p1;
    PLOTAP = line(PL1AP(:,1),PL1AP(:,2));
    hold on
    PLOTA2P = line(p1(:,1),p1(:,2));
    
    
    axis([-0.25 1.25 0 10])
    
        
    elseif gd.MENU == 5   %%%%%%% BETA �LABE 2
        
    CHECK2 = 0;
    delete(PLOT1)
    delete(PLOT2)
    delete(PLOT3)
    delete(PLOT4)
    delete(INLET)
    delete(OUTLET)
    delete(P0)
    delete(P1)
    delete(P2)
    delete(P3)
    delete(P4)
    delete(P5)
    delete(P6)
    delete(P7)
    delete(P8)
    delete(P9)
    delete(PLOTPBA1)
    delete(PLOTPBS1)
    delete(PLOTPBA2)
    delete(PLOTPBS2)
    delete(PLOTBA)
    delete(PLOTBS)
    delete(PLOT2PBA1)
    delete(PLOT2PBS1)
    delete(PLOT2PBA2)
    delete(PLOT2PBS2)
    delete(PLOT2BA)
    delete(PLOT2BS)
    
    delete(PLOTA)
    delete(PLOTA2)
    delete(A1)
    delete(A2)
    delete(A3)
    delete(A4)
    
    delete(PLOTAP)
    delete(PLOTA2P)
    delete(A1P)
    delete(A2P)
    delete(A3P)
    delete(A4P)
    
    delete(PLOTB)
    delete(PLOTB2)
    delete(B1)
    delete(B2)
    delete(B3)
    delete(B4)
    
    delete(PLOTBP)
    delete(PLOTB2P)
    delete(B1P)
    delete(B2P)
    delete(B3P)
    delete(B4P)
    
    delete(PLOT2A)
    delete(PLOT2A2)
    delete(A12)
    delete(A22)
    delete(A32)
    delete(A42)
    
    delete(PLOT2AP)
    delete(PLOT2A2P)
    delete(A12P)
    delete(A22P)
    delete(A32P)
    delete(A42P)
    
    delete(PLOT2B)
    delete(PLOT2B2)
    
    delete(PLOT2BP)
    delete(PLOT2B2P)

    if CHECKB2 == 1
        
        PB12 = getPosition(B12); %%%%%%%%%% RA�Z
        PB22 = getPosition(B22);
        PB32 = getPosition(B32);
        PB42 = getPosition(B42);
        
        delete(B12)
        delete(B22)
        delete(B32)
        delete(B42)
        
        B12 = impoint(gca,0,PB12(2));
        setString(B12,'R1');
        B22 = impoint(gca,PB22(1),PB22(2));
        setString(B22,'R2');
        B32 = impoint(gca,PB32(1),PB32(2));
        setString(B32,'R3');
        B42 = impoint(gca,1,PB42(2));
        setString(B42,'R4');
        
        PB12P = getPosition(B12P); %%%%%%%%%% PUNTA
        PB22P = getPosition(B22P);
        PB32P = getPosition(B32P);
        PB42P = getPosition(B42P);
        
        delete(B12P)
        delete(B22P)
        delete(B32P)
        delete(B42P)
        
        B12P = impoint(gca,0,PB12P(2));
        setString(B12,'P1');
        B22P = impoint(gca,PB22P(1),PB22P(2));
        setString(B22,'P2');
        B32P = impoint(gca,PB32P(1),PB32P(2));
        setString(B32,'P3');
        B42P = impoint(gca,1,PB42P(2));
        setString(B42,'P4');
        
    elseif CHECKB2 == 0

        if TIPO == 1
        
        B12 = impoint(gca,0,0); %%%%%%%%%% RA�Z
        setString(B12,'R1');
        B22 = impoint(gca,0.25,-7.5);
        setString(B22,'R2');
        B32 = impoint(gca,0.5,-15);
        setString(B32,'R3');
        B42 = impoint(gca,1,-30);
        setString(B42,'R4');
        
        B12P = impoint(gca,0,0); %%%%%%%%%% PUNTA
        setString(B12P,'P1');
        B22P = impoint(gca,0.25,-7.5);
        setString(B22P,'P2');
        B32P = impoint(gca,0.5,-15);
        setString(B32P,'P3');
        B42P = impoint(gca,1,-30);
        setString(B42P,'P4');
        
        else
            
        B12 = impoint(gca,0,30); %%%%%%%%%% RA�Z
        setString(B12,'R1');
        B22 = impoint(gca,0.25,22.5);
        setString(B22,'R2');
        B32 = impoint(gca,0.5,15);
        setString(B32,'R3');
        B42 = impoint(gca,1,0);
        setString(B42,'R4');
        
        B12P = impoint(gca,0,20); %%%%%%%%%% PUNTA
        setString(B12P,'P1');
        B22P = impoint(gca,0.25,15);
        setString(B22P,'P2');
        B32P = impoint(gca,0.5,10);
        setString(B32P,'P3');
        B42P = impoint(gca,1,0);
        setString(B42P,'P4');
        
        end
        
        PB12 = getPosition(B12);%%%%%%%%%%% RA�Z
        PB22 = getPosition(B22);
        PB32 = getPosition(B32);
        PB42 = getPosition(B42);
        
        PB12P = getPosition(B12P);%%%%%%%%%%% PUNTA
        PB22P = getPosition(B22P);
        PB32P = getPosition(B32P);
        PB42P = getPosition(B42P);
        
        CHECKB2 = 1;
        
    end
    
        [p1]= [0 PB12(2);PB22(1) PB22(2);PB32(1) PB32(2);1 PB42(2)]; %%%%%%%%%% RA�Z
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1B2=l*p1;
    PLOT2B = line(PL1B2(:,1),PL1B2(:,2));
    hold on
    PLOT2B2 = line(p1(:,1),p1(:,2));
    
           [p1]= [0 PB12P(2);PB22P(1) PB22P(2);PB32P(1) PB32P(2);1 PB42P(2)]; %%%%%%%%%% PUNTA
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1B2P=l*p1;
    PLOT2BP = line(PL1B2P(:,1),PL1B2P(:,2));
    hold on
    PLOT2B2P = line(p1(:,1),p1(:,2));
    
    axis([-0.25 1.25 -60 90])
    
elseif gd.MENU == 6 %%%%%%%%%%%%%%%%%% ESPESOR �LABE 2
    
    CHECKB2 = 0;
    delete(PLOT1)
    delete(PLOT2)
    delete(PLOT3)
    delete(PLOT4)
    delete(INLET)
    delete(OUTLET)
    delete(P0)
    delete(P1)
    delete(P2)
    delete(P3)
    delete(P4)
    delete(P5)
    delete(P6)
    delete(P7)
    delete(P8)
    delete(P9)
    delete(PLOTPBA1)
    delete(PLOTPBS1)
    delete(PLOTPBA2)
    delete(PLOTPBS2)
    delete(PLOTBA)
    delete(PLOTBS)
    delete(PLOT2PBA1)
    delete(PLOT2PBS1)
    delete(PLOT2PBA2)
    delete(PLOT2PBS2)
    delete(PLOT2BA)
    delete(PLOT2BS)
    
    delete(PLOTA)
    delete(PLOTA2)
    delete(A1)
    delete(A2)
    delete(A3)
    delete(A4)
    
    delete(PLOTAP)
    delete(PLOTA2P)
    delete(A1P)
    delete(A2P)
    delete(A3P)
    delete(A4P)
    
    delete(PLOTB)
    delete(PLOTB2)
    delete(B1)
    delete(B2)
    delete(B3)
    delete(B4)
    
    delete(PLOTBP)
    delete(PLOTB2P)
    delete(B1P)
    delete(B2P)
    delete(B3P)
    delete(B4P)
    
    delete(PLOT2A)
    delete(PLOT2A2)
    
    delete(PLOT2AP)
    delete(PLOT2A2P)
    
    delete(PLOT2B)
    delete(PLOT2B2)
    delete(B12)
    delete(B22)
    delete(B32)
    delete(B42)
    
    delete(PLOT2BP)
    delete(PLOT2B2P)
    delete(B12P)
    delete(B22P)
    delete(B32P)
    delete(B42P)
    
    if CHECK2 == 1
        
        PA12 = getPosition(A12);%%%%%%%%%%%%%%%% RA�Z
        PA22 = getPosition(A22);
        PA32 = getPosition(A32);
        PA42 = getPosition(A42);
        
        delete(A12)
        delete(A22)
        delete(A32)
        delete(A42)
        
        A12 = impoint(gca,PA12(1),PA12(2));
        setString(A12,'R1');
        A22 = impoint(gca,PA22(1),PA22(2));
        setString(A22,'R2');
        A32 = impoint(gca,PA32(1),PA32(2));
        setString(A32,'R3');
        A42 = impoint(gca,PA42(1),PA42(2));
        setString(A42,'R4');
        
        PA12P = getPosition(A12P);%%%%%%%%%%%%%%%% PUNTA
        PA22P = getPosition(A22P);
        PA32P = getPosition(A32P);
        PA42P = getPosition(A42P);
        
        delete(A12P)
        delete(A22P)
        delete(A32P)
        delete(A42P)
        
        A12P = impoint(gca,PA12P(1),PA12P(2));
        setString(A12P,'P1');
        A22P = impoint(gca,PA22P(1),PA22P(2));
        setString(A22P,'P2');
        A32P = impoint(gca,PA32P(1),PA32P(2));
        setString(A32P,'P3');
        A42P = impoint(gca,PA42P(1),PA42(P2));
        setString(A42P,'P4');
        
    elseif CHECK2 == 0
        
        A12 = impoint(gca,0,2); %%%%%%%%%%%%%%%% RA�Z
        setString(A12,'R1');
        A22 = impoint(gca,0.2,6);
        setString(A22,'R2');
        A32 = impoint(gca,0.6,2);
        setString(A32,'R3');
        A42 = impoint(gca,1,1);
        setString(A42,'R4');
        
        PA12 = getPosition(A12);
        PA22 = getPosition(A22);
        PA32 = getPosition(A32);
        PA42 = getPosition(A42);
        
        A12P = impoint(gca,0,2); %%%%%%%%%%%%%%%% PUNTA
        setString(A12P,'P1');
        A22P = impoint(gca,0.2,6);
        setString(A22P,'P2');
        A32P = impoint(gca,0.6,2);
        setString(A32P,'P3');
        A42P = impoint(gca,1,1);
        setString(A42P,'P4');
        
        PA12P = getPosition(A12P);
        PA22P = getPosition(A22P);
        PA32P = getPosition(A32P);
        PA42P = getPosition(A42P);
        
        CHECK = 1;
        
    end
    
    
    [p1]= [PA12(1) PA12(2);PA22(1) PA22(2);PA32(1) PA32(2);PA42(1) PA42(2)]; %%%%%%%%%%%%%%%% RA�Z
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1A2=l*p1;
    PLOT2A= line(PL1A2(:,1),PL1A2(:,2));
    hold on
    PLOT2A2 = line(p1(:,1),p1(:,2));
    
    [p1]= [PA12P(1) PA12P(2);PA22P(1) PA22P(2);PA32P(1) PA32P(2);PA42P(1) PA42P(2)]; %%%%%%%%%%%%%%%% PUNTA
   
    n=4;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1A2P=l*p1;
    PLOT2AP= line(PL1A2P(:,1),PL1A2P(:,2));
    hold on
    PLOT2A2P = line(p1(:,1),p1(:,2));
    
        
    axis([-0.25 1.25 0 10])
    
        
end

function updateAxes2(axesToUse,gd)
axes(axesToUse)

global PL1B PL1A PL1BP PL1AP

global ESPESOR BETA UPBLADE DOWNBLADE CENTRALBLADE PLOTLE PLOTTE
global  TIPO TIPOAlabe1 

delete(UPBLADE)
delete(DOWNBLADE)
% delete(CENTRALBLADE)
delete(PLOTLE)
delete(PLOTTE)

if ESPESOR == 1 & BETA == 1 & TIPOAlabe1 == 1; %%%%%%%%%%%%% RA�Z
    
    %Representaci�n �labe
    
   endPL = length(PL1A);
   endPLB = length(PL1B);
   x = linspace(0,1,100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A(round(endPL/100*i),2);
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
   
   %Arreglo circulos
   
    [p1]= [xarriba(1) centralarriba(1); (x(1)-0.03*y(1)) -0.1 ; xdebajo(1) centraldebajo(1)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1LE=l*p1;

    [p1]= [xarriba(end) centralarriba(end); (x(end)+0.03*y(end)) (central(end)-0.12*y(end)) ; xdebajo(end) centraldebajo(end)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1TE=l*p1;
   
%    CENTRALBLADE = plot(x,central,'b');
%    hold on
   DOWNBLADE = plot(xdebajo,centraldebajo,'b');
   hold on
   UPBLADE = plot(xarriba,centralarriba,'b');
   hold on
   PLOTLE = line(PL1LE(:,1),PL1LE(:,2));
   hold on
   PLOTTE = line(PL1TE(:,1),PL1TE(:,2));
   
   axis([-0.25 1.25 -2 2]) 
   
   elseif ESPESOR == 1 & BETA == 1 & TIPOAlabe1 == 2; %%%%%%%%%%%%% PUNTA
    
    %Representaci�n �labe
    
   endPL = length(PL1AP);
   endPLB = length(PL1BP);
   x = linspace(0,1,100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1AP(round(endPL/100*i),2);
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
   
   %Arreglo circulos
   
    [p1]= [xarriba(1) centralarriba(1); (x(1)-0.03*y(1)) -0.1 ; xdebajo(1) centraldebajo(1)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1LE=l*p1;

    [p1]= [xarriba(end) centralarriba(end); (x(end)+0.03*y(end)) (central(end)-0.12*y(end)) ; xdebajo(end) centraldebajo(end)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1TE=l*p1;
   
%    CENTRALBLADE = plot(x,central,'b');
%    hold on

   DOWNBLADE = plot(xdebajo,centraldebajo,'b');
   hold on
   UPBLADE = plot(xarriba,centralarriba,'b');
   hold on
   PLOTLE = line(PL1LE(:,1),PL1LE(:,2));
   hold on
   PLOTTE = line(PL1TE(:,1),PL1TE(:,2));
   

   
   axis([-0.25 1.25 -2 2]) 
   
   
end
   

function updateAxes3(axesToUse,gd)
axes(axesToUse)


global PL1B2 PL1A2 PL1B2P PL1A2P 
global UPBLADE2 DOWNBLADE2 CENTRALBLADE2 PLOTLE2 PLOTTE2 ESPESOR2 BETA2
global TIPOAlabe2 

delete(UPBLADE2)
delete(DOWNBLADE2)
% delete(CENTRALBLADE2)
delete(PLOTLE2)
delete(PLOTTE2)

if ESPESOR2 == 1 & BETA2 == 1 & TIPOAlabe2 == 1; %%%%%%%%%%%%% RA�Z
    
    %Representaci�n �labe
    
   endPL = length(PL1A2);
   endPLB = length(PL1B2);
   x = linspace(0,1,100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A2(round(endPL/100*i),2);
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
   
   %Arreglo circulos
   
    [p1]= [xarriba(1) centralarriba(1); (x(1)-0.03*y(1)) -0.1 ; xdebajo(1) centraldebajo(1)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1LE=l*p1;

    [p1]= [xarriba(end) centralarriba(end); (x(end)+0.03*y(end)) (central(end)-0.12*y(end)) ; xdebajo(end) centraldebajo(end)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1TE=l*p1;
   
%    CENTRALBLADE2 = plot(x,central,'b');
%    hold on

   DOWNBLADE2 = plot(xdebajo,centraldebajo,'b');
   hold on
   UPBLADE2 = plot(xarriba,centralarriba,'b');
   hold on
   PLOTLE2 = line(PL1LE(:,1),PL1LE(:,2));
   hold on
   PLOTTE2 = line(PL1TE(:,1),PL1TE(:,2));
   axis([-0.25 1.25 -2 2]) 
   
elseif ESPESOR2 == 1 & BETA2 == 1 & TIPOAlabe2 == 2; %%%%%%%%%%%%% PUNTA
    
   %Representaci�n �labe
    
   endPL = length(PL1A2P);
   endPLB = length(PL1B2P);
   x = linspace(0,1,100);
   deltax = x(2)-x(1);
   degtorand = 3.14166/180;
   
   %Ley de espesor
   for i=1:100
       y(i)=PL1A2P(round(endPL/100*i),2);
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
   
   %Arreglo circulos
   
    [p1]= [xarriba(1) centralarriba(1); (x(1)-0.03*y(1)) -0.1 ; xdebajo(1) centraldebajo(1)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1LE=l*p1;

    [p1]= [xarriba(end) centralarriba(end); (x(end)+0.03*y(end)) (central(end)-0.12*y(end)) ; xdebajo(end) centraldebajo(end)];
   
    n=3;
    n1=n-1;
    for    i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];

    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation
    end
    
    PL1TE=l*p1;
   
%    CENTRALBLADE2 = plot(x,central,'b');
%    hold on

   DOWNBLADE2 = plot(xdebajo,centraldebajo,'b');
   hold on
   UPBLADE2 = plot(xarriba,centralarriba,'b');
   hold on
   PLOTLE2 = line(PL1LE(:,1),PL1LE(:,2));
   hold on
   PLOTTE2 = line(PL1TE(:,1),PL1TE(:,2));
   axis([-0.25 1.25 -2 2]) 
end


% --- Outputs from this function are returned to the command line.
function varargout = Bezier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global MENUANT ESPESOR BETA MENUANT2 ESPESOR2 BETA2

if get(hObject,'Value') == 2
    MENUANT = 1;
elseif get(hObject,'Value') == 3
    MENUANT = 2;
    MENUANT2 = 2;
    BETA = 1;
elseif get(hObject,'Value') == 4
    MENUANT = 3;
    MENUANT2 = 3;
    ESPESOR = 1;
elseif get(hObject,'Value') == 5
    MENUANT = 4;
    MENUANT2 = 4;
    BETA2 = 1;
elseif get(hObject,'Value') == 6
    MENUANT = 5;
    MENUANT2 = 5;
    ESPESOR2 = 1;
end

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p1x_Callback(hObject, eventdata, handles)
% hObject    handle to p1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1x as text
%        str2double(get(hObject,'String')) returns contents of p1x as a double
global P1 PP1

delete(P1)
axes(handles.axes1)
P1 = impoint(gca,str2double(get(hObject,'String')),PP1(2));
setString(P1,'P1');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p1x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p2x_Callback(hObject, eventdata, handles)
% hObject    handle to p2x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p2x as text
%        str2double(get(hObject,'String')) returns contents of p2x as a double
global P2 PP2

delete(P2)
axes(handles.axes1)
P2 = impoint(gca,str2double(get(hObject,'String')),PP2(2));
setString(P2,'P2');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p2x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p3x_Callback(hObject, eventdata, handles)
% hObject    handle to p3x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p3x as text
%        str2double(get(hObject,'String')) returns contents of p3x as a double
global P3 PP3

delete(P3)
axes(handles.axes1)
P3 = impoint(gca,str2double(get(hObject,'String')),PP3(2));
setString(P3,'P3');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p3x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p3x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p1y_Callback(hObject, eventdata, handles)
% hObject    handle to p1y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1y as text
%        str2double(get(hObject,'String')) returns contents of p1y as a double
global P1 PP1

delete(P1)
axes(handles.axes1)
P1 = impoint(gca,PP1(1),str2double(get(hObject,'String')));
setString(P1,'P1');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p1y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p3y_Callback(hObject, eventdata, handles)
% hObject    handle to p3y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p3y as text
%        str2double(get(hObject,'String')) returns contents of p3y as a double
global P3 PP3

delete(P3)
axes(handles.axes1)
P3 = impoint(gca,PP3(1),str2double(get(hObject,'String')));
setString(P3,'P3');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p3y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p3y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p44y_Callback(hObject, eventdata, handles)
% hObject    handle to p44y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p44y as text
%        str2double(get(hObject,'String')) returns contents of p44y as a double
global P4 PP4

delete(P4)
axes(handles.axes1)
P4 = impoint(gca,PP4(1),str2double(get(hObject,'String')));
setString(P4,'P4');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p44y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p44y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p2y_Callback(hObject, eventdata, handles)
% hObject    handle to p2y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p2y as text
%        str2double(get(hObject,'String')) returns contents of p2y as a double
global P2 PP2

delete(P2)
axes(handles.axes1)
P2 = impoint(gca,PP2(1),str2double(get(hObject,'String')));
setString(P2,'P2');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p2y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p44x_Callback(hObject, eventdata, handles)
% hObject    handle to p44x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p44x as text
%        str2double(get(hObject,'String')) returns contents of p44x as a double
global P4 PP4

delete(P4)
axes(handles.axes1)
P4 = impoint(gca,str2double(get(hObject,'String')),PP4(2));
setString(P4,'P4');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p44x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p44x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p5y_Callback(hObject, eventdata, handles)
% hObject    handle to p5y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p5y as text
%        str2double(get(hObject,'String')) returns contents of p5y as a double
global P5 PP5

delete(P5)
axes(handles.axes1)
P5 = impoint(gca,PP5(1),str2double(get(hObject,'String')));
setString(P5,'P5');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p5y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p5y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p5x_Callback(hObject, eventdata, handles)
% hObject    handle to p5x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p5x as text
%        str2double(get(hObject,'String')) returns contents of p5x as a double
global P5 PP5

delete(P5)
axes(handles.axes1)
P5 = impoint(gca,str2double(get(hObject,'String')),PP5(2));
setString(P5,'P5');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p5x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p5x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p6y_Callback(hObject, eventdata, handles)
% hObject    handle to p6y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p6y as text
%        str2double(get(hObject,'String')) returns contents of p6y as a double
global P6 PP6

delete(P6)
axes(handles.axes1)
P6 = impoint(gca,PP6(1),str2double(get(hObject,'String')));
setString(P6,'P6');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p6y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p6y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p6x_Callback(hObject, eventdata, handles)
% hObject    handle to p6x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p6x as text
%        str2double(get(hObject,'String')) returns contents of p6x as a double
global P6 PP6

delete(P6)
axes(handles.axes1)
P6 = impoint(gca,str2double(get(hObject,'String')),PP6(2));
setString(P6,'P6');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p6x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p6x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p7y_Callback(hObject, eventdata, handles)
% hObject    handle to p7y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p7y as text
%        str2double(get(hObject,'String')) returns contents of p7y as a double
global P7 PP7

delete(P7)
axes(handles.axes1)
P7 = impoint(gca,PP7(1),str2double(get(hObject,'String')));
setString(P7,'P7');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p7y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p7y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p7x_Callback(hObject, eventdata, handles)
% hObject    handle to p7x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p7x as text
%        str2double(get(hObject,'String')) returns contents of p7x as a double
global P7 PP7

delete(P7)
axes(handles.axes1)
P7 = impoint(gca,str2double(get(hObject,'String')),PP7(2));
setString(P7,'P7');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p7x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p7x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p0x_Callback(hObject, eventdata, handles)
% hObject    handle to p0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0x as text
%        str2double(get(hObject,'String')) returns contents of p0x as a double
global P0 PP0

delete(P0)
axes(handles.axes1)
P0 = impoint(gca,str2double(get(hObject,'String')),PP0(2));
setString(P0,'P0');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p0y_Callback(hObject, eventdata, handles)
% hObject    handle to p0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0y as text
%        str2double(get(hObject,'String')) returns contents of p0y as a double
global P0 PP0

delete(P0)
axes(handles.axes1)
P0 = impoint(gca,PP0(1),str2double(get(hObject,'String')));
setString(P0,'P0');

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function p0y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSET.
function pushbuttonSET_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MENUANT

if MENUANT == 1
    MENUANT = 0;
end

gatherAndUpdate(handles)


% --- Executes on slider movement.
function sliderBA2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function sliderBA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderBA1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function sliderBA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderBS2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function sliderBS2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderBS1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function sliderBS1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
Main_RADIAL()
return



function leanangle_Callback(hObject, eventdata, handles)
% hObject    handle to leanangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of leanangle as text
%        str2double(get(hObject,'String')) returns contents of leanangle as a double
global leanangle 

leanangle = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function leanangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leanangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2BA2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2BA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function slider2BA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2BA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2BA1_Callback(hObject, eventdata, handles)
% hObject    handle to slider2BA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function slider2BA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2BA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2BS2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2BS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function slider2BS2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2BS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2BS1_Callback(hObject, eventdata, handles)
% hObject    handle to slider2BS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function slider2BS1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2BS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function leanangle2_Callback(hObject, eventdata, handles)
% hObject    handle to leanangle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of leanangle2 as text
%        str2double(get(hObject,'String')) returns contents of leanangle2 as a double
global leanangle2

leanangle2 = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function leanangle2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leanangle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOmega_Callback(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOmega as text
%        str2double(get(hObject,'String')) returns contents of editOmega as a double
global PerdidaOMEGA

PerdidaOMEGA = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editOmega_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCt_Callback(hObject, eventdata, handles)
% hObject    handle to editCt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCt as text
%        str2double(get(hObject,'String')) returns contents of editCt as a double
global PerdidaCT

PerdidaCT = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editCt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCh_Callback(hObject, eventdata, handles)
% hObject    handle to editCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCh as text
%        str2double(get(hObject,'String')) returns contents of editCh as a double
global PerdidaCH

PerdidaCH = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editT0_Callback(hObject, eventdata, handles)
% hObject    handle to editT0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editT0 as text
%        str2double(get(hObject,'String')) returns contents of editT0 as a double
global temperatura0

temperatura0 = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editT0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editT0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editP0_Callback(hObject, eventdata, handles)
% hObject    handle to editP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editP0 as text
%        str2double(get(hObject,'String')) returns contents of editP0 as a double
global presion0

presion0 = str2double(get(hObject,'String'))*100000;

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editP0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGasto_Callback(hObject, eventdata, handles)
% hObject    handle to editGasto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGasto as text
%        str2double(get(hObject,'String')) returns contents of editGasto as a double

global gastoM

gastoM = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function editGasto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGasto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double
global omega

omega = str2double(get(hObject,'String'));

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuTIPO.
function popupmenuTIPO_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuTIPO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuTIPO contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuTIPO
global TIPO

if get(hObject,'Value') == 1
    TIPO = 1;
else
    TIPO = 2;
end

gatherAndUpdate(handles)
    


% --- Executes during object creation, after setting all properties.
function popupmenuTIPO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuTIPO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
global TIPOAlabe1

if get(hObject,'Value') == 1
    TIPOAlabe1 = 1;
else
    TIPOAlabe1 = 2;
end

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6
global TIPOAlabe2

if get(hObject,'Value') == 1
    TIPOAlabe2 = 1;
else
    TIPOAlabe2 = 2;
end

gatherAndUpdate(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Axes3()
