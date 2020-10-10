% compare numerical and analytic solutions for the homogeneous model
% reproducing Figure 2 of the related manuscript

% clear;clc

FM=60;dt=0.0003;Tc=0.06;
ts=-Tc/2:dt:Tc/2;
wavelet=(1-2*(pi*FM*(ts-3*dt)).^2).*exp(-(pi*FM*(ts-3*dt)).^2);

%wavelet=diff(wavelet)./diff(ts);%derivative of the source wavelet
nr=15;
nt=602;


for testno=1:2


if testno==1
    MT=sqrt(1/2)*[0 1 0;1 0 0;0 0 0];%DC1
    load homo_iso_dc1_2n;%homo_iso_dc2_2n;%
    recvxw1=Sx;
else
    MT=sqrt(1/2)*[0 0 0;0 0 -1;0 -1 0];%DC2
    load homo_iso_dc2_2n;%homo_iso_dc2_2n;%
    recvxw1=Sx;
end


%% cosine direction
ss=[100 100 200];
for i=1:nr
    rx(i)=0;
    ry(i)=0;
    rz(i)=(i-1)*10;
    r(i)=sqrt((rx(i)-ss(1))^2+(ry(i)-ss(2))^2+(rz(i)-ss(3))^2);
    cosx(i)=(rx(i)-ss(1))/r(i);
    cosy(i)=(ry(i)-ss(2))/r(i);
    cosz(i)=(rz(i)-ss(3))/r(i);
    tp(i)=r(i)/3000;
end
    
cosxyz=zeros(nr,3);
cosxyz(:,1)=cosx';
cosxyz(:,2)=cosy';
cosxyz(:,3)=cosz';

ttp=round((tp/dt));
tts=round((tp/dt*1.67));

%% analytical solutions
% The detailed equations of far-field displacement of P- and S-waves
% in a homogeneous isotropic elastic medium
% eqs. (8)(9) of the related manuscript
d=2.5;
vp=3000;
vs=vp/1.67;
Bp=zeros(1,nr);
Bs=zeros(1,nr);
Ap=zeros(1,nr);
As=zeros(1,nr);
pamp=zeros(1,nr);
samp=zeros(1,nr);
for i=1:nr
    Ap(i)=4*pi*r(i)*d*vp.^3;
    As(i)=4*pi*r(i)*d*vs.^3;
    for ii=1:3
        for jj=1:3
            Bp(i)=Bp(i)+cosx(i)*cosxyz(i,ii)*cosxyz(i,jj)*MT(ii,jj);%cos
            Bs(i)=Bs(i)+(KroDel(sym(1),sym(ii))-cosx(i)*cosxyz(i,ii))*cosxyz(i,jj)*MT(ii,jj);%kronecker and cos
        end
    end
    pamp(i)=Bp(i)/Ap(i);
    samp(i)=Bs(i)/As(i);
            
end

wavep=zeros(nr,nt);
waves=zeros(nr,nt);
for i=1:nr
    for j=1:nt
        if j==ttp(i)
            wavep(i,j)=pamp(i);
        end 
        if j==tts(i)
            waves(i,j)=samp(i);
        end 
    end
    
end

for i=1:nr  
    pwave2(i,:)=conv(wavelet,wavep(i,:));
    swave2(i,:)=conv(wavelet,waves(i,:));
    
    
end

TS=0:dt:0.2403;
for i=1:nr
    pwave3(i,:)=diff(pwave2(i,:))./diff(TS);
    swave3(i,:)=diff(swave2(i,:))./diff(TS);
end


figure
h1=subplot (221);
wiggle(pwave3'+swave3','b','b');


%% numerical solution
hold on;
wiggle(recvxw1(:,1:4:4*nr-3),'r','r');
hold on;

set(gca,'fontsize',16);
set(gca,'fontname','arial');
xlabel('Receiver Number ');
set(gca,'ytick',[0 200 400 600 800],'yticklabel',dt*[0 200 400 600 800]);
ylabel('Time (s)');
title('(a)','position',[7.5,860]);%(c)


pswave=pwave3'+swave3';
pswavez=pswave;



nn1=1;nn2=8;nn3=15;
rec=recvxw1;
h2=subplot(222);
trace1=pswave(:,nn1)/max(pswave(:,nn1));
trace2=rec(:,4*nn1-3)/max(rec(:,4*nn1-3));
rmsErr=sqrt(1/length(trace1)*sum((trace1-trace2).^2))
set(gca,'fontname','arial');
plot(trace1,'b-','linewidth',2.5);
hold on;
plot(trace2,'r--','linewidth',2.5);
text(0.05, 0.08,['RMS Error = ',num2str(roundn(rmsErr,-4))],'fontsize',16);
set(gca,'fontname','arial','fontsize',16);
set(gca,'xtick',[0 200 400 600 800],'xticklabel',[]);

axis([0 800 -1 1])
title('(b)');%(d)


h3=subplot(223);
trace1=pswave(:,nn2)/max(pswave(:,nn2));
trace2=rec(:,4*nn2-3)/max(rec(:,4*nn2-3));
rmsErr=sqrt(1/length(trace1)*sum((trace1-trace2).^2))
plot(trace1,'b-','linewidth',2.5);
hold on;
plot(trace2,'r--','linewidth',2.5)
text(0.05, 0.08,['RMS Error = ',num2str(roundn(rmsErr,-4))],'fontsize',16);
set(gca,'fontname','arial','fontsize',16)
set(gca,'xtick',[0 200 400 600 800],'xticklabel',[]);
legend('Analytical','FD method');

%xlabel('T (s)');
ylabel('Normalized Amplitude')
axis([0 800 -1 1])

h4=subplot(224);
trace1=pswave(:,nn3)/max(pswave(:,nn3));
trace2=rec(:,4*nn3-3)/max(rec(:,4*nn3-3));
rmsErr=sqrt(1/length(trace1)*sum((trace1-trace2).^2))
plot(trace1,'b-','linewidth',2.5);
hold on;
plot(trace2,'r--','linewidth',2.5)
text(0.05, 0.08,['RMS Error = ',num2str(roundn(rmsErr,-4))],'fontsize',16);
set(gca,'fontname','arial','fontsize',16)
set(gca,'xtick',[0 200 400 600 800],'xticklabel',dt*[0 200 400 600 800]);
xlabel('T (s)');%ylabel('Normalized Amplitude')
axis([0 800 -1 1])
set(gcf,'pos',[100 100 1200 650])
set(h1,'pos',[0.08 0.1 0.4 0.8]);
set(h2,'pos',[0.575 0.7 0.4 0.2])
set(h3,'pos',[0.575 0.4 0.4 0.2])
set(h4,'pos',[0.575 0.1 0.4 0.2])





end

%% Kronecker Delta function
function output = KroDel(input1,input2)
if input1 == input2
    output = 1;
else
    output = 0;
end
end

