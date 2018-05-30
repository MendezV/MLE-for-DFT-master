function logp=TC_logli(params , MaxPop_F, MaxPop_M, Nbins,Tframes, hist)

%%IN
%%-params: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins, corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop_F:  maximum observed packing in the system for females
%%-MaxPop_M:  maximum observed packing in the system for males
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%calculates the value of the likelihood function for a given set of
%%parameters and data for our specific model

%%OUT
%%-logp: value of the likelihood function for the parameters and data from
%%the analytic formula

%parameters that will be used in the calculation
MaxPop=(MaxPop_F+1)*(MaxPop_M+1); %% index that separates f's and v's in the params vector
V_F=params(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
V_M=params(MaxPop+Nbins+1:end); %%extracting vexation for males
f=params(1:MaxPop); %%extracting frustration

%%gauge fix
f(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
f(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
f(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix

%%index settup for flattened (row major) joint probability distribution
N_Fpre=((1:(MaxPop_F+1))-1)'; %vector with possible female occupation numbers in the system
N_Mpre=((1:(MaxPop_M+1))-1)'; %vector with possible male occupation numbers in the system
N_F=reshape(ones(MaxPop_M+1,1)*N_Fpre',MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
N_M=reshape(N_Mpre*ones(1,MaxPop_F+1),MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution


%%elements that will be needed for the calculation of the negative log of
%%the probability
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
NexpAv_F=hist*N_F; %experimental average number of female flies at each bin
NexpAV_M=hist*N_M; %experimental average number of male flies at each bin
fmeanexp=hist*f; %experimental average frustration at each bin
%loggammaav=mean(log(gamma(counts+1)),2); %%not necessary for the minimization, could do it but sums to
%a huge number 
logp=Tframes*sum(log(z)+V_F.*NexpAv_F+V_M.*NexpAV_M+fmeanexp);

end