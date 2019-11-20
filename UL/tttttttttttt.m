clc
clear
M=25;
K=3;
direc_theta=[-40 10 60];
c=3e8;
f=2.4e9;
lamda=c/f;
m=0.3;
d=lamda*m;
N=1;
r=M*d/(2*pi);
phase=floor(2*pi*r/lamda);
%构造信号和噪声
fs=37.02e6;
n=(1:N);
snr=30;
P=sqrt(10.^(snr/10));
s=exp(1j*2*pi*f*n/fs);
fai1=exp(1j*2*pi*rand(1,N));
fai2=exp(1j*2*pi*rand(1,N));
fai3=exp(1j*2*pi*rand(1,N));
s1=(s.*fai1)';
s2=(s.*fai2)';
s3=(s.*fai3)';
s=[s1,s2,s3]';
nr=randn(M,N);
ni=randn(M,N);
n=nr+1j*ni;
%计算协方差矩阵
atheta1=P*exp(-1j*phase*cos(direc_theta(1)*pi/180-2*pi*[0:M-1]'/M));
atheta2=P*exp(-1j*phase*cos(direc_theta(2)*pi/180-2*pi*[0:M-1]'/M));
atheta3=P*exp(-1j*phase*cos(direc_theta(3)*pi/180-2*pi*[0:M-1]'/M));
atheta=[atheta1,atheta2,atheta3];
x=atheta*s+n;
Rx=x*x'/N;
%构造变换矩阵T
phase=floor(phase);
F=exp(1j*2*pi*[-phase:1:phase]'*[0:1:M-1]/M);
Ji=(1j.^[[-phase:1:phase]']).*besselj([[-phase:1:phase]'],-phase);
J=diag(Ji);
T=inv(J)*F/M;

y=T*x;
Ry=T*Rx*T';
%求空间谱函数


%单次快拍UCA平滑
diag1=(y(phase+1:end)');
diag2=fliplr((y(1:phase+1)'));
Rfb1=toeplitz(diag2,diag1);


%  
Rfb=smooth_ss(Ry,4,K);
%求空间谱函数
[u,d,v] = svd(Rfb1);
Un1 = u(:,K+1:end);
GG1=Un1*Un1';%单次快拍uca
[u1,d1,v1] = svd(Rfb);
Un = u1(:,K+1:end);
GG=Un*Un';%传统uca


direc_search=[-90:0.1:90];
MM=2*phase+1;
for nn=1:length(direc_search)
        a=exp(1j*(-phase:1:phase)'*direc_search(nn)*pi/180);
        b1=a(1:phase+1,:);
        b=a(1:4,:);
        music(nn)= 1/(b'*GG*b);
        music1(nn)=1/(b1'*GG1*b1);
       
end
doa_after=abs(music);
doa_after1=abs(music1);
% doa_after=20*log(abs(music)/max(abs(music)));
 
plot(direc_search ,doa_after1);
title('MODE SPACE-SMOOTH')
xlabel('Angle/(\circ)')
ylabel('Spectrum/dB')
% hold on;
% plot(direc_search ,doa_after);
legend('q');
xlim([-90,90]);
est_direc=sort( find_peak_theta(doa_after1,direc_search,3));
err_direc = cal_theta_match_error_ul(direc_theta,est_direc);



