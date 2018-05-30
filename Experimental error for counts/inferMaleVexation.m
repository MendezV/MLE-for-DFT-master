function [V_M, CovMatV_M]=inferMaleVexation(counts_M)

Nbins=size(counts_M,1); %total number of bins
MaxPop_M=max(max(counts_M)); %maximum observed packing in the system for Males
delta=zeros(MaxPop_M,Nbins); %%allocating the matrix of kronecker delta functions

for n=0:MaxPop_M
    delta(n+1,:)=mean( counts_M'==n ); %%basically computing the joint histogram of counts but flattened
end

hist=delta';
V_M=log(hist(:,1)./hist(:,2));
V_Mexerror=psi(1,hist(:,2))+psi(1,hist(:,1));
V_Mblock=diag(V_Mexerror);

CovMatV_M=[zeros(size(V_Mblock)),zeros(size(V_Mblock));zeros(size(V_Mblock)),V_Mblock];

end