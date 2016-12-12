function landauDamping
np=30000;
l=2*2*pi;
ng=500;
dx=l/ng;
qm=-1;
q=-4.1888e-04;
m=q/qm;
rb=-q*np/l;
gx=0:dx:(ng-1)*dx;
dt=0.01;
n=1;
k=2*pi*n/l;
[x,v]=initializePhaseSpace(np,l);
x0=x;
rhog=accumulateCharge(x,q,ng,rb,dx);
eg=Efield(rhog,l,dx);
ex=interpolateE(eg,x,dx);
v=velAdj(v,ex,dt,q,m);
x=pertX(x,0.1,l,n);
tmax=15;
xv=0:dt:tmax;
i=1;
histPeEnergy = [];
histKeEnergy = [];
histTeEnergy = [];
histSpectrum = [];
time=0:dt:tmax;
for t=time
    v=incVel(v,q,ex,dt,m);
    x=incPos(x,v,dt,l);
    rhog=accumulateCharge(x,q,ng,rb,dx);
    eg=Efield(rhog,l,dx);
    y = fft(eg)'/length(eg);
    histSpectrum = [histSpectrum 2*abs(y(1:ng/2+1))];
    ex=interpolateE(eg,x,dx);
    xv(i)=v(50);
    i=i+1;
    ke =0.5*m*sum(v.*v);
    pe = 0.5*sum(eg.^2)*l;
    te   =  ke + pe;
    histPeEnergy = [histPeEnergy pe];
    histKeEnergy = [histKeEnergy ke];
    histTeEnergy = [histTeEnergy te];
end
% [~,loc]=findpeaks(xv);
% subplot(311),plot((0:dt:tmax),histSpectrum(2,:)),xlabel('time'),ylabel('amptitude of electric field'),grid on;
% subplot(312),plot((0:dt:tmax),histPeEnergy),xlabel('time'),ylabel('kinetic energy'),grid on;
% subplot(313),plot((0:dt:tmax),histKeEnergy),xlabel('time'),ylabel('Field energy'),grid on;
% fprintf('omega: %d k: %d \n',2*pi/((loc(2)-loc(1))*dt),k);
figure
hold on;
plot(time,histPeEnergy,'r');
% plot(time,histKeEnergy,'b');
figure
hist(v);
figure
% hold on;
semilogy(time,histSpectrum(2,:));
hold on;
vt=1;
WP=(np/l)*q*q/m;
dblen=vt/WP;
gamma=sqrt(pi/8)*(WP/(k*dblen)^3)*exp(-((1/(2*(dblen*k)^2))+1.5));
semilogy(time,histSpectrum(2,1)*exp(-gamma*time),'r');

function [x]=pertX(x,x1,l,n)
x=x+(x1*cos(2*pi*(n/l)*x));

function [v]=incVel(v,q,e,dt,m)
v=v+((q/m)*e*dt);

function [x]=incPos(x,v,dt,l)
x=rem(rem(x+(v*dt),l)+l,l);

function[x,v] = initializePhaseSpace(np,l)
x=(l/(2*np)):l/np:l;
% v=zeros(1,np);
v=randn(1,length(x));     

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
% phik(1)=0;
kappa=(k.*sin(k*dx))./(k*dx);
ek=(-1i*kappa).*phik;
% ek(1)=0;
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
