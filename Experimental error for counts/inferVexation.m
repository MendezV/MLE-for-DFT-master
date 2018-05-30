function [V, CovMatV]=inferVexation(counts)

Nbins=size(counts,1); %total number of bins
MaxPop=max(max(counts)); %maximum observed packing in the system for Males
delta=zeros(MaxPop,Nbins); %%allocating the matrix of kronecker delta functions

for n=0:MaxPop
    delta(n+1,:)=mean( counts'==n ); %%basically computing the joint histogram of counts but flattened
end

hist=delta';
V=log(hist(:,1)./hist(:,2));
Vexerror=psi(1,hist(:,2))+psi(1,hist(:,1));
CovMatV=diag(Vexerror);


end