function PK(rho,IPts)
clc
%usage: PK('rho',[IPts])
tic
sym t;
n0=1; %initial condition, n(t=0)
%NBSR Constants
A=6.98E-04; %average neutron generation time
NGr=6; %number of neutron groups
B=[0.00022,0.00111,0.00107,0.00301,0.00092,0.00032]; %delayed neutron fraction of each group
Beta=sum(B); %Total delayed neutron fraction
La=[0.0125,0.0318,0.109,0.317,1.35,8.64]; %decay constants
vfn=1.954752539E+10; %average fission neutron velocity (200 MeV)
CS=585E-24; %microscopic cross section of U-235
Er=3.215568209E-11; %average recoverable energy per U-235 fission
NV=9.814578723E+23; %atomic number density of U-235 in the core's volume

q=rho=='t';
if max(q)>0
    sym(rho);
else
    str2double(rho);
end
Step=0.01; %step size
ti=(0:Step:max(IPts)); %time vector for the simulation
C0=((B./(A.*La)).*n0).'; %initial delayed neutron precursor densities for all the groups
y0=[n0;C0]; %initial y column vector

[t,y]=ode15s(@(t,y) [((double(subs(rho,t))-Beta)/A) La;
    (B(1)/A) -La(1) 0 0 0 0 0;
    (B(2)/A) 0 -La(2) 0 0 0 0;
    (B(3)/A) 0 0 -La(3) 0 0 0;
    (B(4)/A) 0 0 0 -La(4) 0 0;
    (B(5)/A) 0 0 0 0 -La(5) 0;
    (B(6)/A) 0 0 0 0 0 -La(6)]*y,ti,y0); %solving the stiff ODE
n=y(:,1);
P=n*vfn*CS*Er*NV;

dor=zeros(length(IPts),NGr+2);
Power=zeros(length(IPts),2);
for i=1:length(IPts)
    dire= t==IPts(i);
    dor(i,2:end)=y(dire,:);
    Power(i,2)=P(dire);
end
dor(:,1)=IPts';
Power(:,1)=IPts';
disp('Neutron Densities are below')
disp(dor)
disp('power change is below')
disp(Power)

save('PK_densities.txt','dor','-ascii')
save('PK_power.txt','Power','-ascii')
toc %ends the timer for CPU time
end