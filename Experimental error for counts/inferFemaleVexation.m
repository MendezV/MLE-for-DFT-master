function [V_F, CovMatV_F]=inferFemaleVexation(counts_F)

Nbins=size(counts_F,1); %total number of bins
MaxPop_F=max(max(counts_F)); %maximum observed packing in the system for Males
delta=zeros(MaxPop_F,Nbins); %%allocating the matrix of kronecker delta functions

for n=0:MaxPop_F
    delta(n+1,:)=mean( counts_F'==n ); %%basically computing the joint histogram of counts but flattened
end

hist=delta';
V_F=log(hist(:,1)./hist(:,2));
V_Fexerror=psi(1,hist(:,2))+psi(1,hist(:,1));
V_Fblock=diag(V_Fexerror);

CovMatV_F=[V_Fblock,zeros(size(V_Fblock));zeros(size(V_Fblock)),zeros(size(V_Fblock))];

end