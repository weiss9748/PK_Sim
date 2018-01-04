clear; clf; close all; clc
tic
input= fopen('inpa.csv');
inp=textscan(input,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',',');

co1=inp{1,1}; co2=inp{1,2}; co3=inp{1,3};
co4=inp{1,4}; co5=inp{1,5}; co6=inp{1,6};
co7=inp{1,7}; co8=inp{1,8}; co9=inp{1,9};
co10=inp{1,10}; co11=inp{1,11}; co12=inp{1,12};
co13=inp{1,13};
co=[co1,co2,co3,co4,co5,co6,co7,co8,co9,co10,co11,co12,co13];

%constants
gam=1/2; a21=2; a31=48/25; c21=-8.0; c31=372/25; c32=12/5; c41=-112/125;
c42=-54/125; c43=-2/5; b1=19/9; b2=1/2; b3=25/108; b4=125/108;
e1=17/54; e2=7/36; e3=0; e4=125/108;
c1=1/2; c2=-3/2; c3=121/50; c4=29/250;
a2=1;a3=3/5; tiny=1E-30;

Title=cell2mat(co1(1));
Test_Name=cell2mat(co2(3));
n0=str2double(cell2mat(co2(4)));
A=str2double(cell2mat(co2(5)));
rxnum=str2double(cell2mat(co2(7)));
rxstep=str2double(co(8,2:rxnum+1));
rho=str2double(co(9,2:rxnum+1));
drho=str2double(co(10,2:rxnum+1));
NGr=str2double(cell2mat(co2(12)));
B=str2double(co(13,2:NGr+1));
Beta=sum(B);
La=str2double(co(15,2:NGr+1));
Step_Size=cell2mat(co2(17));
NIntr=str2double(cell2mat(co2(19))); %number of interest points
IPts=str2double(co(20,2:NIntr+1)); %interest pts

dfdy=[((rho-Beta)/A) La; (B(1)/A) -La;
    (B(2)/A) -La;
    (B(3)/A) -La;
    (B(4)/A) -La;
    (B(5)/A) -La;
    (B(6)/A) -La];
dfdt0=[((drho*n0)/A);0;0;0;0;0;0];
C0=((B./(A.*La)).*n0).';
y0=[n0;C0];
f0=dfdy*y0;
I=eye(7);
h=0.01;
cont=1;
t=0;
g=[1:cont];
pq=animatedline('Linewidth',0.5)

if Step_Size=='Y'
    while t<max(IPts)
        [a,b,c]=lu(((1./(gam.*h)*I))-dfdy);
        g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
        g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a);
        g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);
        g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c);
        eps=1E-6;
        err=(e1*g1)+(e2*g2)+(e3*g3)+(e4*g4);
        yscale=abs(y0)+abs(h.*f0)+tiny;
        errmax=max(abs(err./yscale));
        if errmax>eps
            while errmax>eps
                [a,b,c]=lu(((1./(gam.*h)*I))-dfdy);
                g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
                g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a);
                g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);
                g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c);
                eps=1E-6;
                err=(e1*g1)+(e2*g2)+(e3*g3)+(e4*g4);
                yscale=abs(y0)+abs(h.*f0)+tiny;
                errmax=max(abs(err./yscale));
                h=max([(0.9*h*(errmax^(-1/3))), 0, 5*h]);
                cont=cont+1;
                t=t+h;
            end
        end
       for i=length(t)
            if errmax>0.1296
                for g=1:cont
                    h=0.9*h*(errmax^-0.25);
                    [a,b,c]=lu(((1./(gam.*h)*I))-dfdy);
                    g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
                    g2=a*((h.*c2.*dfdt0)+((c21*g1)/h))%/norm(a);
                    g3=b*((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))%/norm(b);
                    g4=c*((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)%/norm(c);
                    y=[y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4))]
                    t=t+h;
                    hold on
                    addpoints(pq,t(g),y(g));
                    scatter(t(g),y(g));
                    drawnow;
                    pause(0.01);
                    xlabel('t/s')
                    ylabel('n(t)')
                end
                
            else
                for g=1:cont
                    h=1.5*h;
                    [a,b,c]=lu(((1./(gam.*h)*I))-dfdy);
                    g1=(f0+(h.*c1.*dfdt0))/norm(((1./(gam.*h)*I))-dfdy);
                    g2=((h.*c2.*dfdt0)+((c21*g1)/h))/norm(a);
                    g3=((h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h))/norm(b);
                    g4=((h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3)))/h)/norm(c);
                    y=[y0+((b1*g1)+(b2*g2)+(b3*g3)+(b4*g4))]
                    t=t+h
                    n=y(1);
                    hold on
                    addpoints(pq,t,n);
                    scatter(t,n);
                    drawnow;
                    pause(0.01);
                    xlabel('t/s')
                    ylabel('n(t)')
                end
            end
        end
    end
elseif Step_Size=='N'
    h=str2double(cell2mat(co2(18)));
    g1=(f0+(h*c1*dfdt0));
    g2=(h.*c2.*dfdt0)+((c21*g1)/h);
    g3=(h.*c3.*dfdt0)+(((c31*g1)+(c32*g2))/h);
    g4=(h.*c4.*dfdt0)+(((c41*g1)+(c42*g2)+(c43*g3))/h);
else
    disp('Automatic Step Size should be indicated as Y or N only!')
end

fprintf('CPU Time: %f \n',toc)
