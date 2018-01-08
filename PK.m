function PK
%clear; clf; close all; clc
diary('PK_Output.txt')
disp('--------------------------------------------------------------------------------------------------------------')
disp('This function was developed to solve the reactor point kinetics equation using the Rosenbrock 4th order method')
disp('--------------------------------------------------------------------------------------------------------------')
file=input('Input File Name (with file extension): ','s'); %imports input file as 'file'
tic %starts timer to calculate CPU time
inpu= fopen(file); %saves input file content into
inp=textscan(inpu,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',','); %scans the input file

%seperates the input file into multiple cells based on columns that are
%delimited using commas from the input comma seperated file (csv)
co1=inp{1,1}; co2=inp{1,2}; co3=inp{1,3};
co4=inp{1,4}; co5=inp{1,5}; co6=inp{1,6};
co7=inp{1,7}; co8=inp{1,8}; co9=inp{1,9};
co10=inp{1,10}; co11=inp{1,11}; co12=inp{1,12};
co13=inp{1,13};
co=[co1,co2,co3,co4,co5,co6,co7,co8,co9,co10,co11,co12,co13];

%Shampine parameters
gam=1/2; a21=2; a31=48/25; c21=-8.0; c31=372/25; c32=12/5; c41=-112/125;
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
rho=str2double(co(9,2:rxnum+1)); %reactivity driving function
drho=str2double(co(10,2:rxnum+1)); %1st derivative of the reactivity func.
NGr=str2double(cell2mat(co2(12))); %number of neutron groups
B=str2double(co(13,2:NGr+1)); %delayed neutron fraction of each group
Beta=sum(B); %Total delayed neutron fraction
La=str2double(co(15,2:NGr+1)); %decay constants
Step_Size=cell2mat(co2(17)); %auto step size response (Y or N)
NIntr=str2double(cell2mat(co2(19))); %number of interest points
IPts=str2double(co(20,2:NIntr+1)); %interest pts

dfdy=[((rho-Beta)/A) La; (B(1)/A) -La;
    (B(2)/A) -La;
    (B(3)/A) -La;
    (B(4)/A) -La;
    (B(5)/A) -La;
    (B(6)/A) -La];
dfdt0=[((drho*n0)/A);0;0;0;0;0;0]; %df/dy
C0=((B./(A.*La)).*n0).'; %initial delayed neutron precursor densities for all the groups
y0=[n0;C0]; %initial y column vector
f0=dfdy*y0; %initial f
I=eye(NGr+1); %identity matrix
h=0.01; %initial tentative stepsize
cont=1; %counter for the amount of times the loop runs for automatic stepsize
t=0; %initial time
g=[1:cont]; %used for the for loop later

%plots for the neutron densities can be displayed in subplots (uncommented)
%or in individual figures (uncomment 

%figure for n(t) and C(t)
FIGS=figure('Name','Neutron & Delayed Neutron Precursor Density','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
%comment line directly above for individual, seperate plots
%figure for the neutron density (n) graph
%FIG1=figure('Name','Neutron Density','NumberTitle','off')
%uncomment line directly above to enable seperate plots
subplot(1,2,1)
pq=animatedline('Color','k','Linewidth',0.5);
title('Neutron Density')
xlabel('t/s')
ylabel('n(t)')
%movegui(FIG1,'west');
%uncomment line directly above for seperate plots

%figure for the delayed precursor neutron density
%FIG2=figure('Name','Delayed Neutron Precursor Density','NumberTitle','off')
%uncomment line directly above to enable seperate plots
subplot(1,2,2)
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
%movegui(FIG2,'east');
%uncomment line directly above for seperate plots

%step size (h) determination, g and y column vectors calculations start
if Step_Size=='Y' %if automatic stepsize is set to Y, it is on
    while t<max(IPts) %while the time the calculation is on is less than the maximum of the interest points
        [a,b,c]=lu(((1./(gam.*h)*I))-dfdy); %LU decompisition of LHS of eqns 3 in the paper, based on h=0.01
        g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy); %g1 initial calculation
        g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a); %g2 initial calculation
        g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);%g3 initial calculation
        g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c); %g4 initial calculation
        eps=1E-6; %epsilon, predefined error control criterion
        err=(e1*g1)+(e2*g2)+(e3*g3)+(e4*g4); %initial truncation error calculations
        yscale=abs(y0)+abs(h.*f0)+tiny; %initial truncation error scale
        errmax=max(abs(err./yscale)); %initial maximum error of the y vector
        if errmax>eps %if the initial maximum error of the y vector is greater than epsilon
            while errmax>eps %repeat the following steps until errmax is less than epsilon
                [a,b,c]=lu(((1./(gam.*h)*I))-dfdy); %LU decompisition based on recent h
                g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy); %new g1
                g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a); %new g2
                g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b); %new g3
                g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c); %new g4
                eps=1E-6;
                err=(e1*g1)+(e2*g2)+(e3*g3)+(e4*g4); %new truncation error
                yscale=abs(y0)+abs(h.*f0)+tiny; %new truncation error scale
                errmax=max(abs(err./yscale)); %new maximum error
                h=max([(0.9*h*(errmax^(-1/3))), 0, 5*h]); %new h
                cont=cont+1; %count of how many times this loop was repeated
                t=t+h; %time at which the y value was calculated, based on the h value
            end
        end %errmax is now less than epsilon
        for i=t %repeat the following steps for the length of the time value
            if errmax>0.1296 %if the errmax is greater than 0.1296 do following steps, otherwise, skip to the else statement
                for g=1:max(IPts) %repeat as many times as the maximum of the interest points
                    h=0.9*h*(errmax^-0.25);  %calculating h_next
                    [a,b,c]=lu(((1./(gam.*h)*I))-dfdy); %LU decompisition again but with new h
                    g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy); %new g1
                    g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a); %new g2
                    g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b); %new g3
                    g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c); %new g4
                    y=[y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4))] %new y
                    t=t+h; %time recorded at the calculated y
                    n=y(1); %neutron density at calculated time, n(t)
                    C1=y(2); %delayed neutron precursor density of 1st group
                    C2=y(3); %... of 2nd group
                    C3=y(4); %... of 3rd group
                    C4=y(5); %4th
                    C5=y(6); %5th
                    C6=y(7); %6th
                    hold on %plotting of the values on the 2 independent figures starts
                    addpoints(pq,t,n);
                    scatter(t,n);
                    drawnow;
                    pause(0.01);
                    addpoints(pqk,t,C1);
                    scatter(t,C1);
                    addpoints(pqr,t,C2);
                    scatter(t,C2);
                    addpoints(pqb,t,C3);
                    scatter(t,C3);
                    addpoints(pqg,t,C4);
                    scatter(t,C4);
                    addpoints(pqy,t,C5);
                    scatter(t,C5);
                    addpoints(pqm,t,C6);
                    scatter(t,C6);
                    drawnow;
                    pause(0.01);
                    hold off
                end
            else
                for g=1:max(IPts) %same condition as previous for loop
                    h=1.5*h; %new h, h_next
                    [a,b,c]=lu(((1./(gam.*h)*I))-dfdy); %same steps as previous for loop
                    g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
                    g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a);
                    g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);
                    g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c);
                    y=[y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4))]
                    t=t+h
                    n=y(1);
                    C1=y(2);
                    C2=y(3);
                    C3=y(4);
                    C4=y(5);
                    C5=y(6);
                    C6=y(7);
                    hold on
                    addpoints(pq,t,n);
                    scatter(t,n);
                    drawnow;
                    pause(0.01);
                    addpoints(pqk,t,C1);
                    scatter(t,C1);
                    addpoints(pqr,t,C2);
                    scatter(t,C2);
                    addpoints(pqb,t,C3);
                    scatter(t,C3);
                    addpoints(pqg,t,C4);
                    scatter(t,C4);
                    addpoints(pqy,t,C5);
                    scatter(t,C5);
                    addpoints(pqm,t,C6);
                    scatter(t,C6);
                    drawnow;
                    pause(0.01);
                    hold off
                end
            end
        end
    end
elseif Step_Size=='N'
    h=str2double(cell2mat(co2(18)));
    for g=1:max(IPts) %same condition as previous for loop
        [a,b,c]=lu(((1./(gam.*h)*I))-dfdy); %same steps as previous for loop
        g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
        g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a);
        g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);
        g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c);
        y=[y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4))]
        t=t+h
        n=y(1);
        C1=y(2);
        C2=y(3);
        C3=y(4);
        C4=y(5);
        C5=y(6);
        C6=y(7);
        hold on
        addpoints(pq,t,n);
        scatter(t,n);
        drawnow;
        pause(0.01);
        addpoints(pqk,t,C1);
        scatter(t,C1);
        addpoints(pqr,t,C2);
        scatter(t,C2);
        addpoints(pqb,t,C3);
        scatter(t,C3);
        addpoints(pqg,t,C4);
        scatter(t,C4);
        addpoints(pqy,t,C5);
        scatter(t,C5);
        addpoints(pqm,t,C6);
        scatter(t,C6);
        drawnow;
        pause(0.01);
        hold off
    end
else %if the input file did not specify whether or not automatic step size is activated
    disp('Automatic Step Size should be indicated as Y or N only!') %display this message
end %end step size, g and y vector calculations program

comptime=toc; %tabulate how long it took to do the computation (CPU time)
fprintf('CPU time: %f \n',toc) %print CPU time
end