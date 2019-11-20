function [GG] = smooth_smd(R_in,p,K)
[M,~]=size(R_in); 

q=M+1-p;
J=fliplr(eye(p,p));
Rf=zeros(p,p);
Rb=zeros(p,p);

        for i=1:q
        Rf_this=R_in(i:i+p-1,i:i+p-1);
        Rb_this=J*conj(Rf_this)*J;
        Rf=Rf+Rf_this;
        Rb=Rb+Rb_this;
        end 
        
R=[Rf,Rb];
[U,~,~]=svd(R);
Un=U(:,K+1:p);
GG=Un*Un';   

end

