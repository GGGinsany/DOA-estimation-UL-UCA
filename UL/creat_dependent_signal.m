function X=creat_dependent_signal(A,N,snr)
[~,K]=size(A);
freq=randperm(K+10,K);
S=2*exp(1i*(freq'*(1:N)));
X=A*S;
X=X+awgn(X,snr);
end

