%==============�����ģ����Բ��ͨ���л������ռ�music�㷨
%==============
clear;
close all;
%*****************������������***************************%
 
M=10;

% K_range=1:20;
% for i=1:length(K_range)
%     K=K_range(i)
N=40 ;
K=1;
d2pi=pi/180; 
% pitch_theta=[50 50 70];
% direc_theta=[50 20 100];              %Բ����������˵���Էֱ�360��������
pitch_theta=50*ones(1 ,K);
direc_theta=sort(((2*rand (1,K))-1)*180); 
scan=[-180:1:180]*d2pi;  
thetas = [0:1:90]*d2pi;
lamda=0.3;
k0=2*pi/lamda;                          %����
radius= lamda;                          %����r=��
phase_mode=floor( k0*radius);           %Բ���ܼ��������ģʽ��
MM=2*phase_mode+1;
A=zeros(M,K);
snr=15;
Rm=2*pi*(0:M-1)/M;                      %Բ����Ԫ���ڵĽǶ�
c=3e8;
T=1e-06;
%********************���ɵ���ʸ��*************************%
 
for t=1:K
A(:,t)=exp(1j*k0*(radius*sin(pitch_theta(t)*d2pi)*cos(direc_theta(t)*d2pi*ones(M,1)-Rm')-c*T));
end
 
%*********************�����ռ�ת��***********************%

F = exp(1i*2*pi*(-phase_mode:1:phase_mode)'*(0:1:M-1)/M);
J = (1i.^(-phase_mode:1:phase_mode)').*besselj((-phase_mode:1:phase_mode)',-k0);
Fs=inv(diag(J))*F/M;
% B=Fs*A;  
 
%********************���ɷ����ź�*************************%
 
X=creat_dependent_signal(A,N,snr);
% X=creat_coherent_signal(B,N,snr);  
X=Fs*X;
 Ry = X*X'./N; 
%  Ry=Ry(1:5,1:5);
 Ry=Ry(2:6,2:6);
%*********************music�㷨���***********************%

[u1,d1,v1] = svd(Ry);
Un1 = u1(:,K+1:end);
GG=Un1*Un1'; 
 sp_mode = zeros(length(scan),1);
        for pscan=1:length(scan)
            for tscan = 1:length(thetas)
 
                sigma_temp = sin(thetas(tscan));
                a = exp(1i*k0*  (radius*sigma_temp*cos(Rm'-ones(M,1)*scan(pscan))-c*T)        );
                b = Fs*a;
                b=b(2:6,1);
%                 sp_mode(pscan,tscan)=1/norm(b'*Un1);
                sp_mode(pscan,tscan)=abs((b'*b) /(b'*GG*b) );
            end
        end
 [Z,Y] = meshgrid(thetas/d2pi, scan/d2pi);
 mesh(Z,Y,sp_mode);
% sp_mode = zeros(length(scan),1);
% for  j=1:length(scan)
%     scan=scan(j);
%     v=exp(1i*(-phase_mode:1:phase_mode)*scan);
%     sp_mode(j)=1/norm((v*Rn*GG),2);
%     plot(scan,sp_mode);
% end


p_music_pitch=max( sp_mode);
p_music_direc=max( sp_mode,[],2);
% est_pitch=find_peak_theta(p_music_pitch,thetas/d2pi, K );
est_direc=sort(find_peak_theta(p_music_direc,scan/d2pi,K));
err_direc = cal_theta_match_error(direc_theta,est_direc);
% err_pitch= cal_theta_match_error(pitch_theta,est_pitch);
 return   
    
% end

%*********************root-music�㷨���******************%
         
 
     syms z;
    exp_t = [0:MM-1]';
    pz = z.^exp_t;
    pz_minus = (1/z).^[0:MM-1];
    fz = z.^(MM-1)*pz_minus*GG*pz;
    eq = sym2poly( fz );
    ra=roots(eq);
    rb=ra(abs(ra)<1);                      
 
    [~,I]=sort( abs(rb), 'descend');
 
    w = angle(rb(I(1:K)));
    est_theta=w/pi*180;               %����Ƕ�ֵ
    est_theta=sort(est_theta);
    %disp(est_theta);
   
   return      
 
%*******************��������ͼ**********************%
err = cal_theta_match_error(direc_theta,est_theta);
 
