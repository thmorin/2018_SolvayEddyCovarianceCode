%% Plug flow Model of DO and BOD
% Analytical Solution forward Euler
clear;close;clc;
%% BOD
x=0:0.5:50; %distance and timestamp
dx=0.5; %timestamp
D=nan(size(x)); D(1)=1; %mg L^-1, Deficit
L=nan(size(x)); L(1)=12; %mg L^-1, inital conditions
U=0.3; %m^3 s^-1
H=3; %m
W=10; %m
kd=0.2; %day^-1
ka=3.9*U^0.5/H^1.5; %day^-1
U=0.3*60*60*24; %m^-3 d^-1
U=U/H/W; %m^-3 d^-1
for i=1:length(x)-1
 L(i+1)=L(i)+dx*(-kd/U*L(i));
 D(i+1)=D(i)+dx*( kd/U*L(i)-ka/U*D(i));
end
DO_sat=9;
DO=DO_sat-D;
D_true=kd*L(1)/(ka-kd)*(exp(-kd/U*x)-exp(-ka/U*x))+D(1)*exp(-
ka/U*x);
DO_true=DO_sat-D_true;
L_true=L(1)*exp(-kd/U*x);
DO_err=abs(DO_true-DO);
L_err =abs(L_true -L);
f=figure();
f.Units='Normalized';
f.Position=[0 0 1 1];

%% Plotting
subplot(3,1,1);
h=plotyy(x,L,x,DO);
ylabel(h(1),'L_{approx} (mg L^{-1})');
ylabel(h(2),'DO_{approx} (mg L^{-1})');
xlabel('Distance (m)');
subplot(3,1,2);
h=plotyy(x,L_true,x,DO_true);
ylabel(h(1),'L_{true} (mg L^{-1})');
ylabel(h(2),'DO_{true} (mg L^{-1})');
xlabel('Distance (m)');
subplot(3,1,3);
h=plotyy(x,L_err,x,DO_err);
ylabel(h(1),'L_{error} (mg L^{-1})');
ylabel(h(2),'DO_{error} (mg L^{-1})');
xlabel('Distance (m)');