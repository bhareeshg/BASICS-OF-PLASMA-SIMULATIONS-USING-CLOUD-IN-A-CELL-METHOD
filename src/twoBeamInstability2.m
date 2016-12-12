function coldPlasmaNew2
% hold on;
np=4000;
l=2*pi;
ng=2000;
dx=l/ng;
qc=-0.056;
q=qc*ones(1,np);
q(2:2:np)=-1*q(1:2:np);
m=1;
% rb=-qc*np/l;
rb=0;
gx=0:dx:(ng-1)*dx;
dt=0.01;
n=1;
k=2*pi*n/l;
x=initializeX(np,l);
x0=x;
v=initializeV(np);
rhog=accumulateCharge(x,q,ng,rb,dx);
eg=Efield(rhog,l,dx);
ex=interpolateE(eg,x,dx);
v=velAdj(v,ex,dt,q,m);
x=pertX(x,0.001,l,n);
tmax=40;
xv=0:dt:tmax;
i=1;
for t=0:dt:tmax
%     title('phase space of two beam instability of opposite charges with parameters, np=4000,l=2*pi, ng=2000 ,dt =0.01 and q=-0.056');
    if(t==14)
        subplot(221),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=14'); 
    elseif(t==15)
        subplot(222),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=15');
    elseif(t==16)
        subplot(223),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=16');
    elseif(t==33)
        subplot(224),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=33');
    end
    v=incVel(v,q,ex,dt,m);
    x=incPos(x,v,dt,l);
    rhog=accumulateCharge(x,q,ng,rb,dx);
    eg=Efield(rhog,l,dx);
    ex=interpolateE(eg,x,dx);
    xv(i)=v(50);
    i=i+1;
end
[~,loc]=findpeaks(xv);
figure
plot((0:dt:tmax),xv);
fprintf('omega: %d k: %d \n',2*pi/((loc(2)-loc(1))*dt),k);

function [x]=pertX(x,x1,l,n)
x=x+(x1*cos(2*pi*(n/l)*x));

function [v]=incVel(v,q,e,dt,m)
v=v+(((q./m).*e)*dt);

function [x]=incPos(x,v,dt,l)
x=rem(rem(x+(v*dt),l)+l,l);

function[x]= initializeX(np,l)
% x=(l/(2*np)):l/np:l;
x=zeros(1,np);
x(1:2:np)=(l/(2*np)):2*l/np:l;
x(2:2:np)=(l/(2*np)):2*l/np:l;

function[v] = initializeV(np)
v=zeros(1,np);
v(1:2:np)=-1;
v(2:2:np)=1;
 
function[rhog]=accumulateCharge(x,q,ng,rb,dx)
rhog=rb*ones(1,ng);
np=length(x);
for i=1:1:np
    fl=floor(x(i)/dx);
    rhog(fl+1)=rhog(fl+1)+q(i)*(fl+1-(x(i)/dx));
    if fl+1==ng
    rhog(1)=rhog(1)+q(i)*((x(i)/dx)-fl);
    else
    rhog(fl+2)=rhog(fl+2)+q(i)*((x(i)/dx)-fl);
    end
end
rhog=rhog/dx;

function [e]=Efield(rho,L,dx)
rhok=fft(rho);
N=length(rhok);
k=2*pi*(0:1:N-1)/L;
Ksq=(k.*(sin(k*dx/2)./(k*dx/2))).^2;
phik=rhok./Ksq;
phik(1)=0;
kappa=(k.*sin(k*dx))./(k*dx);
ek=(-1i*kappa).*phik;
ek(1)=0;
e=ifft(ek);
e=real(e);

function [ex]=interpolateE(eg,x,dx)
ex=zeros(1,length(x));
ng=length(eg);
for i=1:1:length(x)
    fl=floor((x(i)/dx));
    lf=eg(fl+1)*(fl+1-(x(i)/dx));
    if fl+1==ng
        gf=eg(1)*((x(i)/dx)-fl);
    else
        gf=eg(fl+2)*((x(i)/dx)-fl);
    end
    ex(i)=gf+lf;
end

function [v]=velAdj(v,ex,dt,q,m)
v=v-((2*q.*ex)/(m*dt));
