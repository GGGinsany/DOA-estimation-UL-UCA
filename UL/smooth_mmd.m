function [GG] = smooth_mmd(R_in,p,K)
[M,~]=size(R_in); 
l=M+1-p;
J=fliplr(eye(p,p));
J1=fliplr(eye(M,M));
RM1=[];
RM2=[];
for i=1:l
    RM_this=[R_in(i:i+p-1,:)];
    RM_b_this=J*conj(RM_this)*J1;
    RM1=[RM1,RM_this];
    RM2=[RM2,RM_b_this];
end
R=[RM1,RM2];
[U,~,~]=svd(R);
 
Un=U(:,K+1:p);
GG=Un*Un';   
end
         

