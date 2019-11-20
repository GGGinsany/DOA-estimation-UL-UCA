 
clc;
clear;
format long %The data show that as long shaping scientific
K=5;
direc_theta=(rand(1,K)-1)*90; %Direction of arrival
 
W=ones(K,1)*4;%Frequency
 
K=length(W); %信源个数K=4
lambda= 0.3;%波长
d=lambda/2;%二分之一波长 天线间距
M_range=20:20;
N_range= 201:10:201;
for row=1:length(M_range)
    M=M_range(row);
    for col=1:length(N_range)
        N=N_range(col);
direc_theta=[-40.1 -20.09 10.36 30.65 50.9]; %Direction of arrival        
A=zeros(M,K); %To creat a matrix with P row and M column
c=3e8;
T=1e-06;
for i=1:K
A(:,i)=exp(-1i*2*pi/lambda*(d*sin(direc_theta(i)/180*pi)-c*T)*[0:M-1]); %求A矩阵
end

S=2*exp(1i*(W*[1:N])); %模拟信号4x200 
X=A*S;
X=X+awgn(X,10);
R=X*X'/N;
% %%%%----------------SMD-------------------%%%%%  
p=M-5;
GG=smooth_smd(R,p,K);   
MUSICP=[];
for n=-90:0.01:90

a=exp(-1i*2*pi/lambda*(d*sin(n/180*pi)-c*T)*[1:p ]'); 
% MUSICP=[MUSICP,1/(a'*G*G'*a)];
MUSICP=[MUSICP,1/(a'*GG*a)];
MUSICP=real(MUSICP);
end
 
figure,plot(-90:0.01:90,MUSICP),axis([-90,90,-50,inf]),title('MUSIC space soothing算法')
grid on;hold on;
est_direc=sort( find_peak_theta(MUSICP,-90:0.01:90,K));
err_direc(row,col) = cal_theta_match_error_ul(direc_theta,est_direc);
    end
end
