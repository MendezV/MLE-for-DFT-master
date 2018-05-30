function [hist, histerr, ExpAv, meanerr]=experimentAnalysis(counts,tau)

%%IN
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe

%%Calculates the inferred probabilities and their averages as calculated by
%%the expectation values. The error is reported noticing that the values of
%%the probabilities follow a Dirichlet distribution.

%%OUT
%%-prob: NbinsxMaxPop matrix that contains the histograms of the occupations for each bin
%%-proberr: square root of the variance on the values of prob according to
%%the Dirichlet distriution
%%-ExpAv: Nbins column vector, Average occupation in each bin
%%-meanerr: Nbins column vector, error on the average occupation in each
%%bin taking into account correlations up to first order. Assuming the
%%values in prob are randomly sampled from a Dirichlet
%%distribution

MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of frames
N=((1:(MaxPop+1))-1)';  %vector with possible occupation numbers in the system

delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end

hist=delta'; %%NbinsxMaxPop matrix, each row correspond to the probability distribution in each bin
histerr=sqrt(hist.*(1-hist)./(Tframes+1)); %%standard error on the parameters of the probability parameters, calculated from the Dirichlet distribution
ExpAv=hist*N; %%Average occupation in each bin
correrr=zeros(Nbins,1); %%the correlation contribution to the error in the mean to first order
for i=1:Nbins
    correrr(i)=-sum(sum(hist(i,:)'*hist(i,:)./(Tframes+1)-diag(hist(i,:)'*hist(i,:))./(Tframes+1),1),2); %%summing covariances calculated from the Dirichlet distribution
meanerr=sqrt((histerr.^2)*(N.^2)+correrr);%%error on the average occupation in each bin assuming perfect counting



end
