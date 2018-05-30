function [ hist, histerr, CovMatHist, probmat, proberr, CovMatProb ] = Histograms( counts , f, V, CovMat,tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
hist=delta';

N=((1:(MaxPop+1))-1)'; %%possible bin occupations for the dataset that we are trying to predict
z=sum(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
probmat=(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1

CovMatHist=zeros(Nbins,MaxPop+1,MaxPop+1);
histerr=sqrt(hist.*(1-hist)./(Tframes+1));

for i=1:Nbins
    CovMatHist(i,:,:)=hist(i,:)'*hist(i,:)./(Tframes+1)-diag(diag(hist(i,:)'*hist(i,:))./(Tframes+1))+diag(histerr(i,:));
end

CovMatProb=zeros(Nbins,MaxPop+1,MaxPop+1);

%%%%%%%

%%%%%%%
proberr=zeros(size(hist));

for i=1:Nbins
    z=sum(exp(-(V(i))*N'-f')./gamma((N+1)')); %%size Nbinsx1 normalization for probaility in each bin for our model
    NsqensAv=sum(((N.^2)').*exp(-(V(i))*N'-f')./gamma((N+1)'))./z; %%ensemble average of n^2 in our model
    NensAv=sum((N').*exp(-(V(i))*N'-f')./gamma((N+1)'))./z; %%ensemble average of our model
    partialV=-(NsqensAv-NensAv.*NensAv);
    HessFF=-(diag(probmat(i,:))-probmat(i,:)'*probmat(i,:)); %%this hessian is the fisher information matrix 
    HessVF=-probmat(i,:).*(N'- NensAv*ones(1,MaxPop+1));
    Jac=[HessFF(2:end,2:end),HessVF(2:end)';HessVF(2:end),partialV];
    RestCovMat=[CovMat(1:MaxPop,1:MaxPop),CovMat(MaxPop+i,1:MaxPop)';CovMat(MaxPop+i,1:MaxPop),CovMat(MaxPop+i,MaxPop+i)];
    
    CovMatProb(i,:,:)=Jac*RestCovMat*Jac';
    
    proberr(i,:)=diag(Jac*RestCovMat*Jac')';
end



end

