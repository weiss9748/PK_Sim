%function PKM
clear; clf;close all; clc

disp('--------------------------------------------------------------------------------------------------------------')
disp('This function was developed to solve the reactor point kinetics equation using the Rosenbrock 4th order method')
disp('--------------------------------------------------------------------------------------------------------------')
file=input('Input File Name (with file extension): ','s'); %imports input file as 'file'
tic %starts timer to calculate CPU time
inpu= fopen(file); %saves input file content intoinp=textscan(inpu,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',','); %scans the input file
inp=textscan(inpu,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',','); %scans the input file
sym t;

%seperates the input file into multiple cells based on columns that are
%delimited using commas from the input comma seperated file (csv)
co1=inp{1,1}; co2=inp{1,2}; co3=inp{1,3};
co4=inp{1,4}; co5=inp{1,5}; co6=inp{1,6};
co7=inp{1,7}; co8=inp{1,8}; co9=inp{1,9};
co10=inp{1,10}; co11=inp{1,11}; co12=inp{1,12};
co13=inp{1,13};
co=[co1,co2,co3,co4,co5,co6,co7,co8,co9,co10,co11,co12,co13];

%extracting data from the input file,
%the input file data is first recognized as cells, and is then converted to
%strings (using cell2mat command), and then, if necessary, to a double percision float
%(using the str2double command).
Title=cell2mat(co1(1)); %title
Test_Name=cell2mat(co2(3)); %test name
n0=str2double(cell2mat(co2(4))); %initial condition, n(t=0)
A=str2double(cell2mat(co2(5))); %average neutron generation time
rxnum=str2double(cell2mat(co2(7))); %number of reactivity steps
rxstep=str2double(co(8,2:rxnum+1)); %reactivity step(s)
rfun=cell2mat(co2(9));
if rfun=='Y'
    rho=sym(co(10,2:rxnum+1)); %reactivity driving function
    drho=diff(rho); %str2double(co(10,2:rxnum+1)); %1st derivative of the reactivity func.
elseif rfun=='N'
    rho=str2double(cell2mat(co(10,2:rxnum+1))); %reactivity driving function
    drho=0;
else
    disp('Is a function of time or not (Line 9 of input file)')
end
NGr=str2double(cell2mat(co2(13))); %number of neutron groups
B=str2double(co(14,2:NGr+1)); %delayed neutron fraction of each group
Beta=sum(B); %Total delayed neutron fraction
La=str2double(co(16,2:NGr+1)); %decay constants
Step_Size=cell2mat(co2(18)); %auto step size response (Y or N)
NIntr=str2double(cell2mat(co2(20))); %number of interest points
IPts=str2double(co(21,2:NIntr+1)); %interest pts

increment=0.01;
ti=(0:increment:max(IPts)); %time vector for the simulation
C0=((B./(A.*La)).*n0).'; %initial delayed neutron precursor densities for all the groups
y0=[n0;C0]; %initial y column vector
tip=zeros(1,length(ti)*max(rxnum));
rhot=zeros(1,length(ti)*max(rxnum));
g=0;
i=1;


if rxnum>1 % the if rxnum part is under development
    while max(tip)<max(IPts)
        g=g+1
        ix = find(IPts>=tip(i),1);
        while tip(i)<IPts(ix);
            rhot(i)=double(subs(rho(g),tip(i)));
            tip(i)=tip(i)+increment;
            i=i+1;
        end
    end
    [t,y]=ode23s(@(t,y) [((rhot(tip==t)-Beta)/A) La;
        (B(1)/A) -La(1) 0 0 0 0 0;
        (B(2)/A) 0 -La(2) 0 0 0 0;
        (B(3)/A) 0 0 -La(3) 0 0 0;
        (B(4)/A) 0 0 0 -La(4) 0 0;
        (B(5)/A) 0 0 0 0 -La(5) 0;
        (B(6)/A) 0 0 0 0 0 -La(6)]*y,ti,y0); %solving the stiff ODE
elseif rxnum<1
    disp('Number of reactivity steps must be greater than 1')
else
    [t,y]=ode23s(@(t,y) [((double(subs(rho,t))-Beta)/A) La;
        (B(1)/A) -La(1) 0 0 0 0 0;
        (B(2)/A) 0 -La(2) 0 0 0 0;
        (B(3)/A) 0 0 -La(3) 0 0 0;
        (B(4)/A) 0 0 0 -La(4) 0 0;
        (B(5)/A) 0 0 0 0 -La(5) 0;
        (B(6)/A) 0 0 0 0 0 -La(6)]*y,ti,y0); %solving the stiff ODE
end

%figure for the neutron density (n) graph
FIG1=figure('Name','Neutron Density','NumberTitle','off');
hold on
title('Neutron Density')
plot(t,y(:,1))
xlabel('t/s')
ylabel('n(t)')
movegui(FIG1,'west');
hold off

%figure for the neutron density (n) graph
FIG2=figure('Name','Neutron Precursor Density','NumberTitle','off');
hold on
title('Neutron Precursor Density')
plot(t,y(:,2),t,y(:,3),t,y(:,4),t,y(:,5),t,y(:,6),t,y(:,7))
xlabel('t/s')
ylabel('C_i(t)')
legend('C_1(t)','C_2(t)','C_3(t)','C_4(t)','C_5(t)','C_6(t)')
movegui(FIG2,'east');
hold off

dor=zeros(length(IPts),NGr+2);
for i=1:length(IPts)
    dire=find(t==IPts(i));
    dor(i,2:end)=y(dire,:);
end
dor(:,1)=IPts';
disp(dor)

toc %ends the timer for CPU time
%end