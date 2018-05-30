function [logratio,D]=modelcomp(params1,params2,CovMat1,CovMat2,counts,tau)
%%IN
%%params1: MaxPop+1+Nbins values for the parameters that minmize the MLE
%%with frustration
%%params2: Nbins values for the parameters that minmize the MLE
%%without frustration
%%CovMat1: MaxPop-1+Nbins square matrix for the parameters in the model
%%without changing the gauge to avoid sigularities in CovMat1
%%CovMat2: Nbins square matrix for the parameters in the naive model
%%counts: NbinsxTframes dataset with the occupations in each bin at each
%%independent time frame

%%computes the logarithm of the ratio of the probability of each model being correct by computing the
%%bayesian evidence as the product of the likelihood function for the MLE
%%optimal parameters and the occam factor that amounts to the size change
%%in the distribution after the data comes under the Gaussian approximation

%%OUT
%%logratio:the logarithm of the ratio of evidences for each model



%%%%%%%%
MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
hist=delta';


loglikelihood1=-(logli(params1 , MaxPop, Nbins, Tframes, hist)) %%loglikelihood for DFT model
loglikelihood2=-(logliNaive(params2 ,  MaxPop, Nbins, Tframes, hist))%%log likelihood for the poisson model


log_prior_1=-0.5*0.01*(params1'*params1)-0.5*(MaxPop+1+Nbins)*log(2*pi*100)%% logarithm of the uniform prior distribution for the optimal parameters for the MLE in the DFT model
log_prior_2=-0.5*0.01*(params2'*params2)-0.5*(Nbins)*log(2*pi*100) %% logarithm of the uniform prior distribution for the optimal parameters for the MLE in the Poisson model

log_post1=0.5*log(det(2*pi*CovMat1)) %% Logartihm of the Occam factor as the difference of the logarith of the width of the prior and posterior distribution(Gaussian approximation) for the DFT model
log_post2=0.5*log(det(2*pi*CovMat2)) %% Logartihm of the Occam factor as the difference of the logarith of the width of the prior and posterior distribution(Gaussian approximation) for the Poisson model

logevidence1=loglikelihood1/log(10)+(log_prior_1-log_post1)/log(10) %%Logarithm of the evidence for the DFT model
logevidence2=loglikelihood2/log(10)+(log_prior_2-log_post2)/log(10) %%Logarithm of the evidence for the Poisson model

logratio=(logevidence1-logevidence2);
D=2*(loglikelihood1-loglikelihood2)

end