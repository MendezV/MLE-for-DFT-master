function logp=logli(params , MaxPop,Nbins,Tframes, hist)

%%IN
%%-params: a vector of size MaxPop+1 + Nbins (MaxPop+1 values for the frustration,the number of flies 
%%can go from zero to the maximum observed packing in all the bins. And
%%Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%calculates the value of the likelihood function for a given set of
%%parameters and data for our specific model

%%OUT
%%-logp: value of the likelihood function for the parameters and data from
%%the analytic formula

%parameters that will be used in the calculation
V=params(MaxPop+2:end); %%extracting vexation
f=params(1:MaxPop+1); %%extracting frustration
f(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
f(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system

%%elements that will be needed for the calculation of the negative log of
%%the probability
z=sum(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%vector that contains the partition function of each bin, size Nbinsx1
NexpAv=hist*N; %experimental average number of flies at each bin
fmeanexp=hist*f; %experimental average frustration at each bin
loggammaav=hist*log(gamma(N+1)); %%not necessary for the minimization, could do it but sums to
%a huge number

logp=Tframes*sum(log(z)+V.*NexpAv+fmeanexp+loggammaav)+sum(sum(gammaln(Tframes*hist+1)))-Nbins*gammaln(Tframes+1);


end