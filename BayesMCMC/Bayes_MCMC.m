function [ fhist, Vhist, fbins, Vbins, MaxLike, fmcmc, fmcmc_mode, Vmcmc, Vmcmc_mode, Likewalk, fmcmcError, VmcmcError, CovMatmcmc] = Bayes_MCMC( counts, tau, Niter )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
tic
MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
histo=delta';


xi=rand([1,MaxPop+1+Nbins]); %%root of the sampling algorithm
xi(1)=0; %%Possible Constraints
xi(2)=0; %%Possible Constraints
walks=zeros([Niter+1,MaxPop+1+Nbins]); %%will contain the samples of the parameters for further processing
lwalk=zeros([Niter+1,1]);

walks(1,:)=xi;
lwalk(1)=loglikelihood(xi' , MaxPop, Nbins, Tframes, histo);

for i=1:Niter

    xprime=normrnd(walks(i,:),0.005);
    xprime(1)=0; %%Possible Constraints
    xprime(2)=0; %%Possible Constraints
    loginit=loglikelihood(walks(i,:)' , MaxPop, Nbins, Tframes, histo);
    logprime=loglikelihood(xprime' , MaxPop, Nbins, Tframes, histo);
    alpha=logprime-loginit;
    if alpha>=0.0
        walks(i+1,:)=xprime;
        lwalk(i+1)=logprime;
    else
        beta=log(rand)/log(10);
        if beta<=alpha
            walks(i+1,:)=xprime;
            lwalk(i+1)=logprime;
        else
            walks(i+1,:)=walks(i,:);
            lwalk(i+1)=lwalk(i);
        end
        
    end
    
   
end

bins=30;

delta=zeros(MaxPop-1,bins); %%allocating the matrix of histograms for the frustration distributions
gamma=zeros(MaxPop-1,bins); %%allocating the matrix of bin edges for the frustration edges
fmcmcError=zeros(MaxPop+1,1);
for n=0:MaxPop
    [ftshist,edges] = hist(walks(200000:end,n+1)',bins);
    fmcmcError(n+1)=std(walks(200000:end,n+1),1);
    %ftshist=hist(walks(:,n+1)',bins);
    delta(n+1,:)=ftshist; 
    gamma(n+1,:)=edges(1:end-1);
end
fhist=delta';
fbins=gamma';


bins=30;
delta=zeros(Nbins,bins); %%allocating the matrix of histograms for the vexation distributions
gamma=zeros(Nbins,bins); %%allocating the matrix of bin edges for the vexation edges
VmcmcError=zeros(Nbins,1);
for n=1:Nbins
    [Vtshist,edges] = hist(walks(200000:end,n+MaxPop+1)',bins);
    VmcmcError(n)=std(walks(200000:end,n+MaxPop+1),1);
    %ftshist=hist(walks(:,n+MaxPop+1)',bins);
    delta(n,:)=Vtshist; 
    gamma(n,:)=edges(1:end-1);
end
Vhist=delta';
Vbins=gamma';

%%Maximum Likelihood
[MaxLike,MaxInd]=max(lwalk);
MaxLikeParams=walks(MaxInd,:)';
Vmcmc=MaxLikeParams(MaxPop+2:end); %%extracting vexation at the minimum
fmcmc=MaxLikeParams(1:MaxPop+1); %%extracting frustration at the minimum
Likewalk=lwalk;


%%%Mode estimation

%%%%%%%%%% Vexation modes
[MaxValsV,MaxIndsV]=max(Vhist);
Vmcmc_mode=zeros([Nbins,1]);
for i=1:Nbins
    Vmcmc_mode(i)=Vbins(MaxIndsV(i),i);
end
%%%%%%%%%%% frustration modes
[MaxValsF,MaxIndsF]=max(fhist);
fmcmc_mode=zeros([MaxPop+1,1]);
for i=0:MaxPop
    fmcmc_mode(i+1)=fbins(MaxIndsF(i+1),i+1);
end


%%%%Covariance matrix

%CovMatmcmc=cov(walks(200000:end,3:end));
CovMatmcmc=cov(walks(100000:end,:));
toc


end

