clc
clear
close all;
%=================参数设置
%=================参数设置
M_range=50:5:100;
m_range=0.8:0.01:2; 
err_direc=zeros(length(M_range),length(m_range));
for row=1:length(M_range)
    for col=1:length(m_range)
M=M_range(row);
m=m_range(col);
                                       % r/lamda
K=3;
n=40;                                           % 角度分辨力
p=4;
direc_theta=[-100:n:(-100+n*(K-1))];        %生成固定间隔的角度

T=1e-06;
c=3e8;
f=2.4e9;
lamda=c/f;
N=600;
r=m*lamda;
k0=2*pi/lamda;
phase=floor(2*pi*r/lamda);
Rm=2*pi*(0:M-1)/M; 
d2pi=pi/180;
snr=40;
MM=2*phase+1;

density=M/(2*pi*r);
%=================构造信号

A=zeros(M,K);
for i=1:K
A(:,i)=exp(-1j*k0*(r*cos(direc_theta(i)*pi/180-Rm')-c*T));
end
x=creat_coherent_signal(A,N,snr);

% x=creat_dependent_signal(A,N,snr);
Rx=x*x'/N;

%=================空间转换

F=exp(1j*2*pi*[-phase:1:phase]'*[0:1:M-1]/M);
J=(1j.^((-phase:1:phase)')).*besselj(((-phase:1:phase)'),-phase);
T=inv(diag(J))*F/M;
Ry=T*Rx*T';
% [u1,d1,v1] = svd(Ry);
% Un1 = u1(:,K+1:end);
% GG=Un1*Un1'; 
%=================空间平滑
% y=T*x;
% GG=smooth_toeplitz (Ry,K);

GG = smooth_smd(Ry,p,K);
%=================music测角
direc_search=(-180:0.1:180);
music=zeros(1,length(direc_search));
for nn=1:length(direc_search)

    
    %=================线阵测角
       a=exp(-1j*(-phase:1:phase)'*direc_search(nn)*pi/180);
       b=a(1:p,:);
       music(nn)=1/(b'*GG*b);
       
    %=================二维测角
%        a=exp(-1j*phase*cos(direc_search(nn)*pi/180-Rm'));
%        b=T*a;
%        music(nn)=1/(b'*GG*b);
end
doa_after=abs(music);
% doa_after=20*log(abs(music)/max(abs(music)));


% =================做图比较
% hold on;
% plot(direc_search,doa_after);
% xlabel('Angle/(\circ)')
% ylabel('Spectrum/dB')
% grid on;

est_direc=sort( find_peak_theta(doa_after,direc_search,K));
err_direc(row,col) = cal_theta_match_error_ul(direc_theta,est_direc);
% scatter(density,err_direc(row,col));hold on;
    end
end

