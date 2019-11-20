function X=creat_coherent_signal(A,N,snr)
[~,K]=size(A);

freq=ones(1,K)*50;
S=2*exp(1i*(freq'*(1:N)));
X=A*S;
X=X+awgn(X,snr);
end

