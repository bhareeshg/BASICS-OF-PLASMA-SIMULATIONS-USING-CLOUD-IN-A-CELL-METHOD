function coldPlasmaNew
% hold on;
np=4000;
l=2*pi;
ng=2000;
dx=l/ng;
q=-0.056;
m=1;
rb=-q*np/l;
gx=0:dx:(ng-1)*dx;
dt=0.2;
n=1;
k=2*pi*n/l;
x=initializeX(np,l);
x0=x;
v0=1;
v=initializeV(np,v0);
rhog=accumulateCharge(x,q,ng,rb,dx);
eg=Efield(rhog,l,dx);
ex=interpolateE(eg,x,dx);
v=velAdj(v,ex,dt,q,m);
x=pertX(x,0.001,l,n);
tmax=80;
xv=0:dt:tmax;
i=1;
time=0:dt:tmax;
histSpectrum=[];
for t=time
    
    if(t==0)
        subplot(321),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=0');
    elseif(t==16)
        subplot(322),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=16'); 
    elseif(t==17)
        subplot(323),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=17');
    elseif(t==18)
        subplot(324),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=18');
    elseif(t==34)
        subplot(325),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=34');
    elseif(t==60)
        subplot(326),hold on,grid on,scatter(x(1:2:np),v(1:2:np),'b'),scatter(x(2:2:np),v(2:2:np),'r'),xlabel('x'),ylabel('v'),title('phase space plot at t=60');
    end
    v=incVel(v,q,ex,dt,m);
    x=incPos(x,v,dt,l);
    rhog=accumulateCharge(x,q,ng,rb,dx);
    eg=Efield(rhog,l,dx);
        y=fft(eg)';
        histSpectrum=[histSpectrum 2*abs(y(1:ng/2+1))];
    ex=interpolateE(eg,x,dx);
    xv(i)=v(50);
    i=i+1;
end
[~,loc]=findpeaks(xv);
figure
plot(time,xv);
fprintf('omega: %d k: %d \n',2*pi/((loc(2)-loc(1))*dt),k);
wp=(np/(2*l))*q*q/m;
omg=abs(sqrt(((k*v0)^2) + (wp^2)-(wp*sqrt((4*(k*v0)^2)+wp^2))));
figure;
semilogy(time,histSpectrum(n,:)*1e17);
hold on;
semilogy(time,100*exp(omg*time),'r');
xlabel('time');
ylabel('growth rate');



function [x]=pertX(x,x1,l,n)
x=x+(x1*cos(2*pi*(n/l)*x));

function [v]=incVel(v,q,e,dt,m)
v=v+((q/m)*e*dt);

function [x]=incPos(x,v,dt,l)
x=rem(rem(x+(v*dt),l)+l,l);

function[x]= initializeX(np,l)
% x=(l/(2*np)):l/np:l;
x=zeros(1,np);
x(1:2:np)=(l/(2*np)):2*l/np:l;
x(2:2:np)=(l/(2*np)):2*l/np:l;


function[v] = initializeV(np,v0)
v=zeros(1,np);
v(1:2:np)=-v0;
v(2:2:np)=v0;
 
function[rhog]=accumulateCharge(x,q,ng,rb,dx)
rhog=rb*ones(1,ng);
np=length(x);
for i=1:1:np
    fl=floor(x(i)/dx);
    rhog(fl+1)=rhog(fl+1)+q*(fl+1-(x(i)/dx));
    if fl+1==ng
    rhog(1)=rhog(1)+q*((x(i)/dx)-fl);
    else
    rhog(fl+2)=rhog(fl+2)+q*((x(i)/dx)-fl);
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
v=v-((2*q*ex)/(m*dt));
