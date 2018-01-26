clear; clf;close all; clc

inpu= fopen('inpa.csv'); %saves input file content into
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

%Shampine parameters
gam=1/2; a21=2; a31=48/25; a32=0.24; c21=-8.0; c31=372/25; c32=12/5; c41=-112/125;
c42=-54/125; c43=-2/5; b1=19/9; b2=1/2; b3=25/108; b4=125/108;
e1=17/54; e2=7/36; e3=0; e4=125/108;
c1=1/2; c2=-3/2; c3=121/50; c4=29/250;
a2=1;a3=3/5; tiny=1E-30;

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
h=zeros(1,991);
h(1)=0.01; %initial tentative stepsize
%ti=(0:h:max(IPts)); %initial time
ti=zeros(1,991); %dummy time vector
i=1; %counter for the amount of times the loop runs for automatic stepsize
I=eye(NGr+1); %identity matrix
eps=1E-6; %epsilon, predefined error control criterion
C0=((B./(A.*La)).*n0).'; %initial delayed neutron precursor densities for all the groups
y0=[n0;C0]; %initial y column vector
yall=zeros(NGr+2,length(ti));
yall(:,1)=[0;y0];
y=ones(NGr+1,1);
dfdy0=[((double(subs(rho,ti(1)))-Beta)/A) La; (B(1)/A) -La;
    (B(2)/A) -La;
    (B(3)/A) -La;
    (B(4)/A) -La;
    (B(5)/A) -La;
    (B(6)/A) -La];
f0=dfdy0*y0; %initial f

%figure for n(t) and C(t)
%FIGS=figure('Name','Neutron & Delayed Neutron Precursor Density','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
%comment line directly above for individual, seperate plots
%figure for the neutron density (n) graph
FIG1=figure('Name','Neutron Density','NumberTitle','off');
%uncomment line directly above to enable seperate plots
%subplot(1,2,1)
pq=animatedline('Color','k','Linewidth',0.5);
title('Neutron Density')
xlabel('t/s')
ylabel('n(t)')
movegui(FIG1,'west');
%uncomment line directly above for seperate plots

%figure for the delayed precursor neutron density
FIG2=figure('Name','Delayed Neutron Precursor Density','NumberTitle','off');
%uncomment line directly above to enable seperate plots
%subplot(1,2,2)
pqk=animatedline('Color','k','Linewidth',0.5);
pqr=animatedline('Color','r','Linewidth',0.5);
pqb=animatedline('Color','b','Linewidth',0.5);
pqg=animatedline('Color','g','Linewidth',0.5);
pqy=animatedline('Color','y','Linewidth',0.5);
pqm=animatedline('Color','m','Linewidth',0.5);
title('Delayed Neutron Precursor Density')
legend('C1','C2','C3','C4','C5','C6')
xlabel('t/s')
ylabel('C_i(t)')
movegui(FIG2,'east');
%uncomment line directly above for seperate plots

while max(ti)<max(IPts)
    i=i+1;
    ti(i)=ti(i-1)+h(i-1);
    dfdy=[((double(subs(rho,ti(i)))-Beta)/A) La; (B(1)/A) -La;
        (B(2)/A) -La;
        (B(3)/A) -La;
        (B(4)/A) -La;
        (B(5)/A) -La;
        (B(6)/A) -La];
    dfdt=[((double(subs(drho,ti(i)))*(yall(2,i)))/A);0;0;0;0;0;0]; %df/dt
    [L,U]=lu(((1./(gam*h(i-1))))-dfdy); %LU decompisition of LHS of eqns 3 in the paper, based on h=0.01
    g1=(f0+(h(i-1).*c1.*dfdt))/norm(L*U); %g1 initial calculation
    g2=((([((double(subs(rho,(a2*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
        (B(2)/A) -La;
        (B(3)/A) -La;
        (B(4)/A) -La;
        (B(5)/A) -La;
        (B(6)/A) -La])*(y0+(a21*g1)))+(h(i-1).*c2.*dfdt)+((c21*g1)/h(i-1)))/norm(L*U); %g2 initial calculation
    g3=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
        (B(2)/A) -La;
        (B(3)/A) -La;
        (B(4)/A) -La;
        (B(5)/A) -La;
        (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c3.*dfdt)+(((c31*g1)+(c32*g2))/h(i-1)))/norm(L*U);%g3 initial calculation
    g4=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
        (B(2)/A) -La;
        (B(3)/A) -La;
        (B(4)/A) -La;
        (B(5)/A) -La;
        (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c4.*dfdt)+(((c41*g1)+(c42*g2)+(c43*g3)))/h(i-1))/norm(L*U); %g4 initial calculation
    err=(e1*g1)+(e2*g2)+(e3*g3)+(e4*g4); %initial truncation error calculations
    yscale=abs(y0)+abs(h(i-1).*f0)+tiny; %initial truncation error scale
    errmax=max(abs(err./yscale)); %initial maximum error of the y vector
    if errmax>eps
        h(i)=max([(0.9*h(i-1)*(errmax^(-1/3))), 0, 5*h(i-1)]); %new h
        ti(i)=ti(i-1)+h(i-1);
        dfdy=[((double(subs(rho,ti(i)))-Beta)/A) La; (B(1)/A) -La;
            (B(2)/A) -La;
            (B(3)/A) -La;
            (B(4)/A) -La;
            (B(5)/A) -La;
            (B(6)/A) -La];
        dfdt=[((double(subs(drho,ti(i)))*(yall(2,i)))/A);0;0;0;0;0;0]; %df/dt
        [L,U]=lu(((1./(gam*h(i-1))))-dfdy); %LU decompisition of LHS of eqns 3 in the paper, based on h=0.01
        g1=(f0+(h(i-1).*c1.*dfdt))/norm(L*U); %g1 initial calculation
        g2=((([((double(subs(rho,(a2*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
            (B(2)/A) -La;
            (B(3)/A) -La;
            (B(4)/A) -La;
            (B(5)/A) -La;
            (B(6)/A) -La])*(y0+(a21*g1)))+(h(i-1).*c2.*dfdt)+((c21*g1)/h(i-1)))/norm(L*U); %g2 initial calculation
        g3=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
            (B(2)/A) -La;
            (B(3)/A) -La;
            (B(4)/A) -La;
            (B(5)/A) -La;
            (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c3.*dfdt)+(((c31*g1)+(c32*g2))/h(i-1)))/norm(L*U);%g3 initial calculation
        g4=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
            (B(2)/A) -La;
            (B(3)/A) -La;
            (B(4)/A) -La;
            (B(5)/A) -La;
            (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c4.*dfdt)+(((c41*g1)+(c42*g2)+(c43*g3)))/h(i-1))/norm(L*U); %g4 initial calculation
    elseif errmax<eps
        if errmax>0.1296
            h(i)=0.9*h(i-1)*(errmax^-0.25);
            ti(i)=ti(i-1)+h(i-1);
            dfdy=[((double(subs(rho,ti(i)))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La];
            dfdt=[((double(subs(drho,ti(i)))*(yall(2,i)))/A);0;0;0;0;0;0]; %df/dt
            [L,U]=lu(((1./(gam*h(i-1))))-dfdy); %LU decompisition of LHS of eqns 3 in the paper, based on h=0.01
            g1=(f0+(h(i-1).*c1.*dfdt))/norm(L*U); %g1 initial calculation
            g2=((([((double(subs(rho,(a2*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a21*g1)))+(h(i-1).*c2.*dfdt)+((c21*g1)/h(i-1)))/norm(L*U); %g2 initial calculation
            g3=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c3.*dfdt)+(((c31*g1)+(c32*g2))/h(i-1)))/norm(L*U);%g3 initial calculation
            g4=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c4.*dfdt)+(((c41*g1)+(c42*g2)+(c43*g3)))/h(i-1))/norm(L*U); %g4 initial calculation
        else
            h(i)=1.5*h(i-1);
            ti(i)=ti(i-1)+h(i-1);
            dfdy=[((double(subs(rho,ti(i)))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La];
            dfdt=[((double(subs(drho,ti(i)))*(yall(2,i)))/A);0;0;0;0;0;0]; %df/dt
            [L,U]=lu(((1./(gam*h(i-1))))-dfdy); %LU decompisition of LHS of eqns 3 in the paper, based on h=0.01
            g1=(f0+(h(i-1).*c1.*dfdt))/norm(L*U); %g1 initial calculation
            g2=((([((double(subs(rho,(a2*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a21*g1)))+(h(i-1).*c2.*dfdt)+((c21*g1)/h(i-1)))/norm(L*U); %g2 initial calculation
            g3=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c3.*dfdt)+(((c31*g1)+(c32*g2))/h(i-1)))/norm(L*U);%g3 initial calculation
            g4=((([((double(subs(rho,(a3*h(i-1))))-Beta)/A) La; (B(1)/A) -La;
                (B(2)/A) -La;
                (B(3)/A) -La;
                (B(4)/A) -La;
                (B(5)/A) -La;
                (B(6)/A) -La])*(y0+(a31*g1)+(a32*g2)))+(h(i-1).*c4.*dfdt)+(((c41*g1)+(c42*g2)+(c43*g3)))/h(i-1))/norm(L*U); %g4 initial calculation
        end
    end
    y=y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4));
    yall(:,i)=[ti(i);y];
    n=y(1);
    C1=y(2);
    C2=y(3);
    C3=y(4);
    C4=y(5);
    C5=y(6);
    C6=y(7);
    hold on
    addpoints(pq,ti(i),n);
    %scatter(ti(i),n);
    drawnow;
    %pause(0.01);
    addpoints(pqk,ti(i),C1);
    %scatter(ti(i),C1);
    addpoints(pqr,ti(i),C2);
    %scatter(ti(i),C2);
    addpoints(pqb,ti(i),C3);
    %scatter(ti(i),C3);
    addpoints(pqg,ti(i),C4);
    %scatter(ti(i),C4);
    addpoints(pqy,ti(i),C5);
    %scatter(ti(i),C5);
    addpoints(pqm,ti(i),C6);
    %scatter(ti(i),C6);
    drawnow;
    %pause(0.01);
    hold off
end