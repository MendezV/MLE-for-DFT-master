function [F_TC, V_F, V_M, CovMat_TC, ferror_TC, V_Ferror, V_Merror, Like] = TC_Bayes(counts_M,counts_F,tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


tic

MaxPop_F=max(max(counts_F)); %maximum observed packing in the system for Females
MaxPop_M=max(max(counts_M)); %maximum observed packing in the system for Males
MaxPop=(MaxPop_M+1)*(MaxPop_F+1); %%including the zeroth fly for both genders
gauge=0;

Nbins=size(counts_M,1); %total number of bins
Tframes=size(counts_M,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop,Nbins); %%allocating the matrix of kronecker delta functions

for m=0:MaxPop_F
for n=0:MaxPop_M
delta(n+1+(MaxPop_M+1)*m,:)=mean(counts_M'==n & counts_F'==m ); %%basically computing the joint histogram of counts but flattened
end
end
histo=delta';



dlmwrite('histog1.dat', [Nbins,MaxPop], 'delimiter',' ');
dlmwrite('histog2.dat', histo, 'delimiter',' ');

system("cat histog1.dat histog2.dat >histog.dat");
system("rm histog1.dat histog2.dat ");


system(" gcc Bayes_post_TC.c -o Bayes_TC.x");
strTframes=int2str(Tframes);
n_iter=2000000;
strn_iter=int2str(n_iter);
strMaxPop_M=int2str(MaxPop_M);
strMaxPop_F=int2str(MaxPop_F);

runn= ["./Bayes_TC.x histog.dat " strTframes " " strn_iter " " strMaxPop_M " " strMaxPop_F];
system(runn);

f = csvread("fmle.dat");
V_M = csvread("V_Mmle.dat");
V_F = csvread("V_Fmle.dat");
Like = csvread("like.dat");
F_TC=reshape(f,[MaxPop_M+1,MaxPop_F+1]);



[stderrors,CovMat_TC]=TC_getCovMat(f,V_F,V_M, MaxPop_F,MaxPop_M,Nbins,Tframes,0); %getting the covariance matrix for the parameters in the model without gauge transformation
ferror_TC=[0;0;stderrors(1:1+MaxPop_M);0;stderrors(2+MaxPop_M:end)]; %%unpacking errors
V_Ferror=stderrors(MaxPop-2:MaxPop+Nbins-3); %%unpacking errors for V_F
V_Merror=stderrors(MaxPop+Nbins-2:end);  %%unpacking errors for V_M

toc

end
