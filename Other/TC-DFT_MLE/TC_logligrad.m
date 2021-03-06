function grad=TC_logligrad(params ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist)

%%IN
%%-params: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins, corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop_F:  maximum observed packing in the system for females
%%-MaxPop_M:  maximum observed packing in the system for males
%%-Nbins: total number of bins
%%-Tframes: number of frames


%%calculates the gradient of the likelihood function for our specific model
%%from the analytical formula at the value of the parameters being passed as input with a give set of
%%data

%%OUT
%%-grad: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins that corresponds to the gradient of the likelihood function

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
NexpAv_M=hist*N_M; %experimental average number of male flies at each bin
NensAv_F=sum((ones(Nbins,1)*N_F').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NensAv_M=sum((ones(Nbins,1)*N_M').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of males in our model
Vgrad_F=Tframes*(NexpAv_F-NensAv_F); %V sector of the gradient is the difference between the observed average and the model average for the females
Vgrad_M=Tframes*(NexpAv_M-NensAv_M); %V sector of the gradient is the difference between the observed average and the model average for the females

%frustration sector of the gradient
probmat=(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')))./(z*ones(1,MaxPop)); %size NbinsxMaxPop+1
fgrad=Tframes*(sum(hist,1)-sum(probmat,1))'; % f sector of the gradient as the difference of the model probability and the hitogram
fgrad(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
fgrad(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
fgrad(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
%complete gradient
grad=[fgrad;Vgrad_F;Vgrad_M];

end
