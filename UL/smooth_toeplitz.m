 function [GG] = smooth_toeplitz(R_in,K)
[MM,~]=size(R_in); 
a_diag=zeros(1,MM-1);
a1=sum(diag(R_in,0))/MM;
    for k=1:MM-1
        a_diag(k)=sum(diag(R_in,-k))/(MM-k);
    end
  al_diag=[a1,a_diag];
  au_diag=[a1,conj(a_diag)];
  RT_n=toeplitz(au_diag,al_diag);
  [u1,~,~] = svd(RT_n );         
  Un1 = u1(:,K+1:end);
  GG=Un1*Un1'; 
end
