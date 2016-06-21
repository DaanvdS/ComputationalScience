%% set variables
clear;
close all;
load('ljdata(1).mat');
N = 100;                                                         %number of data points
subplot(2,2,1);
errorbar(r,potr,noisepotr, 'o');
%Hoi
%% Task 1 A: fit the LJ potential
%Matrixalgebra om vector a te bepalen, via matrix B en vector s
s = potr./noisepotr;
g = [r.^-12 r.^-6];                                              %V=a*r met a=[4eo^-12 -4eo^-6]
B = zeros(size(g));                                              %Preallocating for quicker operation
for k = 1:size(g,2)
    B(:,k) = g(:,k)./noisepotr;
end
A = B'*B;
b = B'*s;
a = A\b;

%Eps en sigma bepalen uit: a=[4eo^-12 -4eo^-6]
eps  = (a(2,1)/-4)^2/(a(1,1)/4);
sigm =((a(1,1)/4)/(a(2,1)/-4))^(1/6);

%Make rfit with N values linspaced in same interval as r
r_fit = linspace(r(1,1),r(length(r),1),N)';
Vr    = (4*eps).*((sigm./r_fit).^12-(sigm./r_fit).^6);

hold on;
subplot(2,2,1)
plot(r_fit,Vr,'b');
grid on;
xlabel('Distance r (m)');
ylabel('Energy E (J)');

 
%% Task 1 B determine r_equi of Vr
% Pen & paper led to:
r_equi  = 8.2875*10^(-11);
Vr_equi = -9.1180e-19; 
hold on;
subplot(2,2,1);
plot(r_equi,Vr_equi,'ro');


%% Task 1 C Harmonic oscillator
% set datapoints around r_equi
rh1   = 0.80e-10;                                                   %NB, moet in de buurt van r_equi, anders is het geen goede benadering meer
rh2   = rh1+2*(r_equi-rh1);                                         %Kies RH2 op dezelfde afstand tot r_equi als RH1
rh    = linspace(rh1,rh2,N)';
potrH = (4*eps).*((sigm./rh).^12-(sigm./rh).^6);
noisepotrH = 1e-100.*ones(N,1);                                     %uncertainty van potH is verwaarloosbaar klein
clear s B g A a;                                                    %clear data van oude fit
g = [ rh.^0 (rh-r_equi).^2 ]; 

%zelfde procedure als Task 1A, fout is verwaarloosbaar klein gemaakt
s = potrH./noisepotrH;
B = zeros(size(g));
for k = 1:size(g,2)
    B(:,k)=g(:,k)./noisepotrH;
end
A = B'*B;
b = B'*s;
a = A\b;

%set rh for better plot
rh1 = 0.75e-10;
rh2 = 0.90e-10;
rh  = linspace(rh1,rh2,N)';
P0  = a(1,1);
ks  = 2*a(2,1);
PR  = P0+ks*(rh-r_equi).^2;

hold on;
subplot(2,2,1);
plot(rh,PR,'r');
ylim([-2.5e-18 2.5e-18]);

%% Task 1D 
F_v = -1*dVdr(r_fit,sigm,eps);
F_h = -1*dPdr(r_fit,r_equi,ks);


subplot(2,2,3);
plot(r_fit,F_v);
hold on;
plot(r_equi,-1*dVdr(r_equi,eps,sigm),'ro');
grid on;
xlabel('Distance r (m)');
ylabel('Force F (N)');

subplot(2,2,4);
plot(r_fit,F_h);
hold on;
plot(r_equi,-1*dPdr(r_equi,r_equi,ks),'ro');
grid on;
xlabel('Distance r (m)');
ylabel('Force F (N)');

%% Task 1E
%Taylorreeks maken voor V(r) rond r_equi. De afgeleide van die reeks geeft
%een constante d^2V/dr^2 voor de (r-r_equi) term. Als je dit vergelijkt met
%ks uit 1D, krijg je ks_taylor=d^2V/dt^2 ge-evalueerd in r_equi.
%NB ks_taylor is ong. ks.
ks_taylor = d2Vdr2(r_equi,sigm,eps); 


clearvars -except eps sigm rh r_equi ks

%% Task 2A: Two hydrogen system
%NB: 0K, one-body problem simplification applied.
m1    = 1.67e-27;
m2    = 1.67e-27;
m_eff = m1*m2/(m1+m2);
N     = 500;
tstart=0;   tfinish=1e-14; y1start=r_equi; y2start=1;

    %solve ODE (ordinary differential equations)
    tspan  = linspace(tstart,tfinish,N);
    ystart = [y1start y2start];
    [t,y]  = ode45(@(t,y)hydrogenODE_V(t,y,m_eff,eps,sigm),tspan,ystart);   
x_V = y(:,1);
v_V = y(:,2);
figure;
a=subplot(2,2,1);
plot(t,x_V);  xlabel('t');    ylabel('x_V'); grid on;
subplot(2,2,2);
plot(t,v_V);  xlabel('t');    ylabel('v_V'); grid on;
    clear y t
    
    %solve ODE (ordinary differential equations)
    tspan  = linspace(tstart,tfinish,N);
    ystart = [y1start y2start];
    [t,y]  = ode45(@(t,y)hydrogenODE_P(t,y,m_eff,ks,r_equi),tspan,ystart);   
x_P = y(:,1);
v_P = y(:,2);
subplot(2,2,3);
plot(t,x_P);  xlabel('t');    ylabel('x_P'); grid on;
subplot(2,2,4);
plot(t,v_P);  xlabel('t');    ylabel('v_P'); grid on;

set(gcf,'NextPlot','add');
axes;
h = title('Oscillaties bij V(r) en P(r)');
set(gca,'Visible','off');
set(h,'Visible','on');
%% Task 2B
