 function [GG] = smooth_toeplitz_1(X,K)
 [M,~]=size(X);
R_1=corrcoef(X.');
a_diag=R_1(1,:);
RT_1=toeplitz(a_diag);  
  
%%%%------Toeplitz�Ľ�����2--------%%%%%  


[V1,~]=eig(RT_1);
Un1=V1(:,1:M-K); %��ȡ�����ӿռ� 
GG=Un1*Un1'; 
 end

