function PK(SA1,SA2,SA3,SA4,RR,dV,dTm,IPts)

%WARNING: The contstants used in this script are specifically made for the NBSR.
%WARNING: The void and temperature coefficients assume start-up (SU) conditions in an LEU core.
% USAGE: PK(SA1,SA2,SA3,SA4,RR,dV,dTm,[IPts])
% Usage of PK function:
% Change the working directory for MATLAB to match the directory that the PK file is in.

% The function can be utilized using the following command:

% PK(SA1,SA2,SA3,SA4,RR,dV,dTm,[IPts])

% Where SA1,SA2,SA3,SA4, and RR are the respectively the positions for shim arms 1 through 4 in degrees, the regulating rod position in inches.
% Where dV and dTm are the change in void volume(Liters) and the change in moderator temperature (Celsius) respectively
% Where IPts is a vector representing the time interest points (in seconds) for the simulation

% Here is a guide of how to use the PK function (the one attached):
% Change MATLAB’s working directory to the folder you placed the attached PK script in.
% On the command line, type PK(SA1,SA2,SA3,SA4,RR,dV,dTm,[IPts]). The first five inputs are the shim arms and regulating rod positions respectively, 
% dV and dTm are the change in void volume (Liters) and the change in moderator temperature (Celsius) respectively, 
% and IPts is a vector of your time interest points (the time points, after the reactivity insertion, where you wish to see the output neutron density and power change).
% Feed the command to MATLAB (ie. click Enter)
% Then check the output files (PK_densities.txt shows the neutron densities at the interest points, and PK_power.txt shows the power change at the interest points)

% For additional information, contact:
% - Jorge Medina (jmedina9416@gmail.com)
% - Abdullah Weiss (weiss9748@gmail.com)

tic
syms t x x1 x2 x3 x4;
n0=1; %initial condition, n(t=0)
%NBSR Constants
A=6.98E-04; %average neutron generation time
NGr=14; %number of neutron groups
B=[0.00022,0.00111,0.00107,0.00301,0.00092,0.00032,0.000203,0.000065,0.0000223,0.0000107,0.0000066,0.0000074,0.000001,0.00000033]; %delayed neutron fraction of each group
Beta=sum(B); %Total delayed neutron fraction
La=[0.0125,0.0318,0.109,0.317,1.35,8.64,0.278,0.0169,0.0049,0.00152,4.27E-4,1.16E-4,4.41E-5,3.65E-6]; %decay constants
vfn=1.954752539E+10; %average fission neutron velocity (200 MeV)
CS=585E-24; %microscopic cross section of U-235
Er=3.215568209E-11; %average recoverable energy per U-235 fission
NV=9.814578723E+23; %atomic number density of U-235 in the core's volume
al_V=-0.039; %LEU Void coefficient, all thimbles voided, SU
al_Tm=-0.0280; %LEU moderator temperature coefficient, SU

D1=diff((-0.000000199*x1^5 + 0.000037466*x1^4 - 0.002360314*x1^3 + 0.055344704*x1^2 - 0.104622548*x1 + 0.064537320),1);
D2=diff((0.000000011*x2^6 - 0.000001943*x2^5 + 0.000134803*x2^4 - 0.004732855*x2^3 + 0.077910025*x2^2 - 0.131259930*x2),1) + 0.055881659;
D3=diff((0.000005906*x3^4 - 0.000767074*x3^3 + 0.028367595*x3^2 - 0.095523922*x3 + 0.092542584),1);
D4=diff((0.000015119*x4^4 - 0.001498778*x4^3 + 0.043485231*x4^2 - 0.105207397*x4 + 0.074580293),1);
DR=diff((0.000068446*x^4 - 0.008057260*x^3 + 0.252661726*x^2 + 0.244962887*x),1);
rho_RR=double(subs(DR,RR));
rho_SA=double(subs(D1,SA1))+double(subs(D2,SA2))+double(subs(D3,SA3))+double(subs(D4,SA4));
rho=rho_RR+rho_SA;
q=rho=='t';
if max(q)>0
    sym(rho);
else
    str2double(rho);
end
Step=0.001; %step size for interest points
ti=(0:Step:max(IPts)); %time vector for the simulation
C0=((B./(A.*La)).*n0).'; %initial delayed neutron precursor densities for all the groups
y0=[n0;C0]; %initial y column vector

%void and moderator temperature feedback conditioning
if dTm==0&&dV==0
    feedback=1;
else
    feedback=((al_Tm*dTm)+(al_V*dV));
end

[t,y]=ode15s(@(t,y) [(((double(subs(rho,t))-Beta)/A)*feedback) La;
    (B(1)/A) -La(1) 0 0 0 0 0 0 0 0 0 0 0 0 0;
    (B(2)/A) 0 -La(2) 0 0 0 0 0 0 0 0 0 0 0 0;
    (B(3)/A) 0 0 -La(3) 0 0 0 0 0 0 0 0 0 0 0;
    (B(4)/A) 0 0 0 -La(4) 0 0 0 0 0 0 0 0 0 0;
    (B(5)/A) 0 0 0 0 -La(5) 0 0 0 0 0 0 0 0 0;
    (B(6)/A) 0 0 0 0 0 -La(6) 0 0 0 0 0 0 0 0;
    (B(7)/A) 0 0 0 0 0 0 -La(7) 0 0 0 0 0 0 0;
    (B(8)/A) 0 0 0 0 0 0 0 -La(8) 0 0 0 0 0 0;
    (B(9)/A) 0 0 0 0 0 0 0 0 -La(9) 0 0 0 0 0;
    (B(10)/A) 0 0 0 0 0 0 0 0 0 -La(10) 0 0 0 0;
    (B(11)/A) 0 0 0 0 0 0 0 0 0 0 -La(11) 0 0 0;
    (B(12)/A) 0 0 0 0 0 0 0 0 0 0 0 -La(12) 0 0;
    (B(13)/A) 0 0 0 0 0 0 0 0 0 0 0 0 -La(13) 0;
    (B(14)/A) 0 0 0 0 0 0 0 0 0 0 0 0 0 -La(14)]*y,ti,y0); %solving the stiff ODE
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
disp('Check the output files!')

% disp('Neutron Densities are below')
% disp(dor)
% disp('power change is below')
% disp(Power)

save('PK_densities.txt','dor','-ascii')
save('PK_power.txt','Power','-ascii')
toc %ends the timer for CPU time
end