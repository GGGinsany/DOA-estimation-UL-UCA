%==============这个test文件是用来测试波束空间圆阵测角的性能
%==============现在拿来做微扰恢复的测试
%==============现在拿来做空间平滑的测试
%==============现在用线阵的方法测角
clear;
close all;
%*****************基本参数设置***************************%
% M=20;
% N=2000;
K=1;
d2pi=pi/180;
pitch_theta=50*ones(1 ,K); 
% pitch_theta=sort(rand(1,K)*90); 
% direc_theta=sort(((2*rand (1,K))-1)*180); 
direc_theta=[70];
scan=[-180:1:180]*d2pi; 
%         thetas = [30:1:70]*d2pi;
lamda=0.3;
k0=2*pi/lamda;                       %波数
radius= lamda;                     %设置r=λ
phase_mode=floor(k0*radius);        %圆阵能激励的最大模式数
MM=2*phase_mode+1;
snr=20;
M_range=30;
N_range=2000;
for ce=1:length(M_range)
    M=M_range(ce);
    for d=1:length(N_range)
                N=N_range(d);

        A=zeros(M,K);
        Rm=2*pi*(0:M-1)/M;                %圆阵阵元所在的角度

        %********************生成导向矢量*************************%

        for t=1:K
                     A(:,t)=exp(1j*k0*radius*sin(pitch_theta(t)*d2pi)*cos(direc_theta(t)*d2pi*ones(M,1)-Rm'));
        end

        %*********************波束空间转换***********************%

        F = exp(-1i*2*pi*(-phase_mode:1:phase_mode)'*[0:1:M-1]/M);
        J = (1i.^[[-phase_mode:1:phase_mode]']).*besselj([[-phase_mode:1:phase_mode]'],-k0);
        Fs=inv(diag(J))*F/M;
        

        %********************生成发送信号*************************%

        X=creat_dependent_signal(A,N,snr);
%         X=creat_coherent_signal(B,N,snr);  
        Rx=X*X';
        
     Ry=Fs*Rx*Fs';
%         p=7;
%         q=MM+1-p;
%         J1=fliplr(eye(p,p));
%         Rf=zeros(p,p);
%         Rb=zeros(p,p);
%         for i=1:q
%         Rf_this=Ry(i:i+p-1,i:i+p-1);
%         Rb_this=J1*conj(Rf_this)*J1;
%         Rf=Rf+Rf_this;
%         Rb=Rb+Rb_this;
%         end
%         Ryy=Rf+Rb;






        %*********************music算法测角***********************%
        [u1,d1,v1] = svd(Ry );         
        Un1 = u1(:,K+1:end);
        GG=Un1*Un1'; 

        sp_mode = zeros(length(scan),1);
        sp_mode_ul = zeros(length(scan),1);
                for dscan=1:length(scan)
%                     for pscan = 1:length(thetas)
                        sigma_temp = sin(50*d2pi);
%                         a = exp(1i*k0*sigma_temp*radius*cos(Rm'-ones(M,1)*scan(dscan)));
%                         b =Fs*a; 
%                         b=b(1:p,1);
                          b=exp(1j*(-phase_mode:1:phase_mode)'*scan(dscan));
                        sp_mode(dscan)=abs(1/(b'*GG*b));
         

                      
                end
                plot( scan/d2pi,sp_mode);
%          [X,Y] = meshgrid(thetas/d2pi,scan /d2pi);
%          mesh(X,Y,sp_mode);
        

%         p_music_pitch=max( sp_mode);
%         p_music_direc=max( sp_mode,[],2);
        p_music_direc= sp_mode;
        % est_pitch=find_peak_theta(p_music_pitch,thetas/d2pi, K );
        est_direc=sort(find_peak_theta(p_music_direc,scan/d2pi,K));
        err_direc(ce,d) = cal_theta_match_error(direc_theta,est_direc);
        % err_pitch= cal_theta_match_error(pitch_theta,est_pitch);
   end
end
 return   
    
% end
 
%*********************root-music算法测角******************%
         
 
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
    est_theta=w/pi*180;               %输出角度值
    est_theta=sort(est_theta);
    %disp(est_theta);
   
   return      
 
%*******************误差估计做图**********************%
err = cal_theta_match_error(direc_theta,est_theta);
 