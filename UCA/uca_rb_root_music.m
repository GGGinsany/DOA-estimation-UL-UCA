%============只在俯仰角为50度的时候测量准确 初步估计是等效线阵距离没求对
close all;
d2pi = pi/180;          % 
M = 50;                  % antenna number
test=1:20 ;
K=2;
for i=1:length(test)

doa = sort((rand(K,1) )*90);
  
    lambda=1;                % wavelength
    L = 2000;                 % snapshot
    


    snr = 40;    
    pn = 10^(-snr/20);  % noise power 
    n = pn*(randn(M,L)+1i*randn(M,L))/sqrt(2); %noise
    s = (randn(K,L)+1i*randn(K,L))/sqrt(2);    %signal
    sigma = sin(50*pi/180);

        r =2*lambda;
        beta = 2*pi*r/lambda;
        A = exp(-1i*beta*sigma*cos(2*pi*[0:M-1]'/M*ones(1,K)-ones(M,1)*doa'*d2pi)); 
        x = A*s; % array input


        K0=floor(beta); % maximum value 
        F = exp(1i*2*pi*[-K0:1:K0]'*[0:1:M-1]/M);
        J = (1i.^[[-K0:1:K0]']).*besselj([[-K0:1:K0]'],-beta);
        y = inv(diag(J))*F/M*x;
        Ry = y*y'/L;
        [u1,d1,v1] = svd(Ry);
 

        %============ music
        Un1 = u1(:,K+1:end);
        GG=Un1*Un1'; 

    k0=2*pi/lambda;  
    phase_mode=ceil(k0*r);  
    MM=2*phase_mode-1;


   

    syms z;
    exp_t = [0:MM-1]';
    pz = z.^exp_t;
    pz_minus = (1/z).^[0:MM-1];
    fz = z.^(MM-1)*pz_minus*GG*pz;
    eq = sym2poly( fz );
    ra=roots(eq);
    rb=ra(abs(ra)<1);                        %abs return 绝对值或是复数幅值
                                             % pick the n roots that are closest to the unit circle
    [~,I]=sort( abs(rb), 'descend');
  
    w = angle(rb(I(1:K)));
    est_theta=sort(w) /d2pi;               %输出角度值
 scatter(doa,est_theta);hold on;
end