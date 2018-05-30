function [ Vhist_N, Vbins_N, Vmcmc_N_mode, MaxLike_N, Vmcmc_N, Likewalk_N, VmcmcError_N, CovMatmcmc_N] = Bayes_MCMC_N( counts, tau, Niter )
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


xi=rand([1,Nbins]);
%xi(1)=0;
%xi(2)=0;
walks=zeros([Niter+1,Nbins]);
lwalk=zeros([Niter+1,1]);

walks(1,:)=xi;
lwalk(1)=loglikelihoodN(xi' , MaxPop, Nbins, Tframes, histo);

for i=1:Niter

    xprime=normrnd(walks(i,:),0.001);
    %xprime(1)=0;
    %xprime(2)=0;
    loginit=loglikelihoodN(walks(i,:)' , MaxPop, Nbins, Tframes, histo);
    logprime=loglikelihoodN(xprime' , MaxPop, Nbins, Tframes, histo);
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

delta=zeros(Nbins,bins); %%allocating the matrix of kronecker delta functions
gamma=zeros(MaxPop-1,bins);
VmcmcError_N=zeros(Nbins,1);
for n=1:Nbins
    [Vtshist,edges] = histcounts(walks(100000:end,n)',bins);
    VmcmcError_N(n)=std(walks(100000:end,n),1);
    %ftshist=hist(walks(:,n+MaxPop+1)',bins);
    delta(n,:)=Vtshist; 
    gamma(n,:)=edges(2:end);
end
Vhist_N=delta';
Vbins_N=gamma';

%%Maximum Likelihood
[MaxLike_N,MaxInd]=max(lwalk);
MaxLikeParams=walks(MaxInd,:)';
Likewalk_N=lwalk;
Vmcmc_N=MaxLikeParams; %%extracting vexation at the minimum

%%Mode
%%%%%%%%%% Vexation modes
[MaxValsV,MaxIndsV]=max(Vhist_N);
Vmcmc_N_mode=zeros([Nbins,1]);
for i=1:Nbins
    Vmcmc_N_mode(i)=mean([Vbins_N(MaxIndsV(i),i),Vbins_N(MaxIndsV(i)-1,i)]);
end



%CovMatmcmc=cov(walks(200000:end,3:end));
CovMatmcmc_N=cov(walks(100000:end,:));
toc


end

