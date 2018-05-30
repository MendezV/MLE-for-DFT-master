function [ chisqDFT, chisqNaive, pvalDFT, pvalNaive ] = Chi_squared( counts, f , V, muDFT, VNaive, muNaive, tau)
%%IN
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe
%%-f: a MaxPop+1(MaxPop is the maximum number of flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)
%%-V: a Nbins(number of bins in the system) sized vector that corresponds
%%-VNaive: a Nbins(number of bins in the system) sized vector that
%%corresponds to the vexation that came out of the MLE for the naive model

%%%%%%%%
MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
histo=delta';

%%Pobability for the DFT model
z=sum(exp(-(V-muDFT)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
probmat=(exp(-(V-muDFT)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1


%%Pobability for the Poisson model
z_P=sum(exp(-(VNaive-muNaive)*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
probmat_P=(exp(-(VNaive-muNaive)*N')./gamma(ones(Nbins,1)*(N+1)'))./(z_P*ones(1,MaxPop+1)); %size NbinsxMaxPop+1


c=log(histo);
c(isinf(-log(histo)))=-60;
deg=ones(size(probmat)).*histo./histo;
deg(isnan(deg))=0;
degree=sum(sum(deg))

%%the logarithm f the probability is normally distributed, the weight in
%%the sum is appropiae


%chisqDFT=sum(sum(deg.*((c-log(probmat)).^2)./(-log(probmat)),2)); %%average chi squared in the system per degree of freedom
%chisqNaive=sum(sum(deg.*((c-log(probmat_P)).^2)./(-log(probmat_P)),2));
%pvalDFT=1-gammainc(chisqDFT/2,(degree-(MaxPop+Nbins))/2); %% the number of degrees of freedom in each bin is taken as the maximum observed number of flies (the zeroth value is constrained due to normalization)
%pvalNaive=1-gammainc(chisqNaive/2,(degree-Nbins)/2);


%%the logarithm f the probability is normally distributed, the weight in
%%the sum is not appropiae or the degrees of freedom?
sigmo=psi(1,histo*Tframes)-psi(1,Tframes);

chisqDFT=sum(sum(deg.*((c-log(probmat)).^2)./sigmo),2); %%average chi squared in the system per degree of freedom
chisqNaive=sum(sum(deg.*((c-log(probmat_P)).^2)./sigmo,2));
pvalDFT=1-gammainc(chisqDFT/2,((MaxPop+1)*Nbins-(MaxPop+Nbins))/2); %% the number of degrees of freedom in each bin is taken as the maximum observed number of flies (the zeroth value is constrained due to normalization)
pvalNaive=1-gammainc(chisqNaive/2,((MaxPop+1)*Nbins-Nbins)/2);

%chisqDFT=sum(Tframes*sum(deg.*((histo-probmat).^2)./probmat,2));
%chisqNaive=sum(Tframes*sum(deg.*((histo-probmat_P).^2)./probmat_P,2));
%pvalDFT=1-gammainc(chisqDFT/2,((MaxPop+1)*Nbins-(MaxPop+Nbins-1))/2); %% the number of degrees of freedom in each bin is taken as the maximum observed number of flies (the zeroth value is constrained due to normalization)
%pvalNaive=1-gammainc(chisqNaive/2,((MaxPop+1)*Nbins-(Nbins))/2);

%%the pvalue is the probability that the residual chi-squared is 
%%greater than the observed chi-squared, if it is too small it is because
%%our residual is too large given the number of degrees of freedom so that
%%it is unlikely that the residuals are gaussian and the errors are random.
%%Thus, we reject the null hypothesis that the theoretical probabilities are
%%the ones that produce the data.
end

