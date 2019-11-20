 function [GG] = smooth_toeplitz_2(X,K)
[M,~]=size(X);
R_1=corrcoef(X.');
a_diag=R_1(1,:);

  
%%%%------Toeplitz�Ľ�����2--------%%%%%  

a_diag1=[a_diag(1),conj(a_diag(2:M))];
RT_2=toeplitz(a_diag1,a_diag);  


[V2,~]=eig(RT_2);
Un2=V2(:,1:M-K); %��ȡ�����ӿռ� 
GG=Un2*Un2';
