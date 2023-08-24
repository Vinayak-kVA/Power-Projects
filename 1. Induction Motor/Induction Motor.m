clc;
clear all;
%% Variable declaration
Vm=220*(sqrt(2))/(sqrt(3));
P=4;
N=1710;
F=60;
rs=0.435;
Xls=0.754;
Xm=26.13;
rr=0.816;
Xlr=0.754;
J=0.089;
we=2*pi*F;
wr=0;
Ls=(Xls+Xm)/we;
Lm=Xm/we;
Lr=(Xlr+Xm)/we;
I=zeros(4,1);
delI=zeros(4,1);
V=zeros(4,1);
L=[Ls 0 Lm 0
   0 Ls 0 Lm
   Lm 0 Lr 0
   0 Lm 0 Lr];
LINV=inv(L);
r=[rs 0 0 0
   0 rs 0 0
   0 0 rr 0
   0 0 0 rr];
G=[0 0 0 0
   0 0 0 0
   0 Lm 0 Lr
  -Lm 0 -Lr 0];
C=(2/3)*[1     -1/2       -1/2
         0   sqrt(3)/2 -sqrt(3)/2];
%% Initial Conditions
 t(1)=0;
 wr(1)=0;
 T(1)=0;
 h=0.0001; %Step size 1msec, (1sec/0.1msec= 10000 iterations)
%% Iterations
for x=1:1:30000
%ABC to DQ transformation    
    Vabc(:,x)=[Vm*sin(we*t(x))
               Vm*sin((we*t(x))-(2*pi/3))
               Vm*sin((we*t(x))+(2*pi/3))];
    Vdq(:,x)=C*Vabc(:,x);
    Vds(x)=Vdq(1,x);
    Vqs(x)=Vdq(2,x);
    V(:,x)=[Vds(1,x)
       Vqs(1,x)
         0
         0];
%RK 4th order weight calculation
K1(:,x)= h*(LINV*(V(:,x)-(r*I(:,x))-(wr(x).*G*I(:,x))));
K2(:,x)= h*(LINV*(V(:,x)-(r*(I(:,x)+(K1(:,x)/2)))-(wr(x).*G*(I(:,x)+(K1(:,x)/2)))));
K3(:,x)= h*(LINV*(V(:,x)-(r*(I(:,x)+(K2(:,x)/2)))-(wr(x).*G*(I(:,x)+(K2(:,x)/2)))));
K4(:,x)= h*(LINV*(V(:,x)-(r*(I(:,x)+(K3(:,x))))-(wr(x).*G*(I(:,x)+(K3(:,x))))));
   
% Updating currents
delI(:,x)=(1/6)*(K1(:,x)+(2*K2(:,x))+(2*K3(:,x))+K4(:,x));
I(:,x+1)=I(:,x)+delI(:,x);
    
%Torque Calculation
T(x+1)=(0.75*P*Lm*((I(3,x)*I(2,x))-(I(4,x)*I(1,x))));

% Speed Calculation
wr(x+1)=wr(x)+((P/2)*(T(x+1)/J)*h);

%time update
t(x+1)=t(x)+h;

end
%Mechanical speed in rpm
wrm=wr*(2/P);
Nrm=(wrm*60)/(2*pi);
figure(1);
plot(t,T)
title('Torque Vs Time');
xlabel('Time (Sec)');
ylabel('Torque (N-m)');
figure(2);
plot(wrm,T)
title('Torque Vs Speed');
xlabel('Speed (Rad/Sec)');
ylabel('Torque (N-m)');
figure(3);
plot(Nrm,T)
title('Torque Vs Speed');
xlabel('Speed (RPM)');
ylabel('Torque (N-m)');
figure(4);
plot (t,wrm);
title('Speed Vs Time');
xlabel('Time (Sec)');
ylabel('Speed (Rad/Sec)');
figure(5);
plot (t,Nrm);
title('Speed Vs Time');
xlabel('Time (Sec)');
ylabel('Speed (RPM)');


