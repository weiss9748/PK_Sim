clear; clc

%y=[]; ysave=[]; dydx=[]; dfdy=[]; dfdx=[]; ystart=[]; yscal=[]; err=[];
x1=0.0; initstep=1E-2; eps=1E-6;

%Parameters
safety=0.9; grow=1.5; pgrow=0.25; shrnk=0.5; pshrnk=1/3;
gam=1/2; a21=2; a31=48/25; c21=-8.0; c31=372/25; c32=12/5; c41=-112/125; 
c42=-54/125; c43=-2/5; b1=19/9; b2=1/2; b3=25/108; b4=125/108;
e1=17/54; e2=7/36; e3=0; e4=125/108;
c1x=1/2; c2x=-3/2; c3x=121/50; c4x=29/250;
a2x=1;a3x=3/5; tiny=1E-30;

%solving linear system of eqns
dydtmp=0;
x=x1 %initial x
if autostep==true
    htry=initstep %initial step
else
    htry=fix_stepsize
end

y=ystart; %initial y
flag_out=0;

piece_rho=1; piece_int=1;
disp('Begin Calculation')

IO_data=load('data');
