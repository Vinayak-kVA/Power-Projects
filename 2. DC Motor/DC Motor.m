clc;
clear all;
%% Variable declaration
Va=230;
Vf=230;
P=4;
Ra=2;
La=0.002;
Rf=115;
Md=0.5;
J=0.05;
wr=0;
Ia=0;
delIa=0;
If=Vf/Rf;
%% Initial Conditions
 t(1)=0;
 wr(1)=0;
 T(1)=0;
 Ia(1)=0;
 h=0.0001; %Step size 1msec, (1sec/0.1msec= 10000 iterations)
%% Iterations
for x=1:1:25000
%RK 4th order weight calculation
K1(x)=(h*(Va-(Ra*Ia(x))-(wr(x)*Md*If))/La);
K2(x)=(h*(Va-(Ra*(Ia(x)+(K1(x)/2)))-(wr(x)*Md*If))/La);
K3(x)=(h*(Va-(Ra*(Ia(x)+(K2(x)/2)))-(wr(x)*Md*If))/La);
K4(x)=(h*(Va-(Ra*(Ia(x)+(K3(x))))-(wr(x)*Md*If))/La);
   
% Updating currents
delIa(x)=(1/6)*(K1(x)+(2*K2(x))+(2*K3(x))+K4(x));
Ia(x+1)=Ia(x)+delIa(x);
    
%Torque Calculation
T(x+1)=0.5*P*Md*If*Ia(x);

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