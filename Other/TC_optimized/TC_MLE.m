function [f_TC, V_F, V_M, CovMat_TC, ferror_TC, V_Ferror, V_Merror, TC_MLE_like]=TC_MLE(rootparams,alpharoot,counts_M,counts_F,gauge,tau)
%,CovMat,ferror,Verror
%%%%%%%%%%
%%IN
%%-rootparams: a vector of size  (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins (MaxPop_F(M)+1 values for the frustration,the number of female (male) flies 
%%can go from zero to the maximum observed packing in all the bins. And Nbins values for the vexation at each bin )
%%it corresponds to the root of the minimization algorithm 
%%-rootalpha: scalar that corresponds to the root that will be used in
%%the linesearch algorithm that is implemented in the minimization 
%%-counts_F: a NbinsxTframes matrix with the number of female flies observed in each
%%bin at each timeframe
%%-counts_M: a NbinsxTframes matrix with the number of male flies observed in each
%%bin at each timeframe
%%-gauge numerical value, if 0 no gauge transformation is applied to the
%%parameters and f(0)=f(1)=0. else the gauge is set so that the average
%%potential is zero and f(1) corresponds to the sum of the previous values
%%for the potential
%%-tau:correlation time corrects for the number of independent time samples




%%finds the minimum of g(V,F)=-logP(V,F|data) by the method of preconditioned conjugate
%%gradients the argument of this function at the minimum gives the maximum
%%likelihood estimates for the parameters in our model. The minimization
%%constraints the values of f(0) and f(1) to be 0 so that the covariance
%%matrix is non-singular and there is a well defined gaussian asymptotic
%%limit for the distribution of estimators

%%OUT
%%-f: a (MaxPop_F+1)(MaxPop_M+1) matrix (MaxPop_F(M) is the maximum number of female (male) flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)
%%-V_F: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for females
%%-V_M: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for males
%%-stderrors: a (MaxPop-1+Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the parameters that we are estimating without the
%%gauge fixed values, which have no uncertainty
%%-CovMat: a ((MaxPop_F+1)(MaxPop_M+1)-3+2*Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators. (if there was no gauge fix)
%%-CovMat:a ((MaxPop_F+1)(MaxPop_M+1)-1+2*Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators after performing the gauge transformation
%%which ammounts to performing a similarity transformation to the covriance matrix
%%also we append to the asymptotic covariance the error and covariances of the parameter f(1,0) and f(0,1) that was previously fixed but now has error. (if there was a gauge fix)



%%%%%%%%

more off;

MaxPop_F=max(max(counts_F)); %maximum observed packing in the system for Females
MaxPop_M=max(max(counts_M)); %maximum observed packing in the system for Males
MaxPop=(MaxPop_M+1)*(MaxPop_F+1); %%including the zeroth fly for both genders

Nbins=size(counts_F,1); %total number of bins
Tframes=size(counts_F,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop,Nbins); %%allocating the matrix of kronecker delta functions

for m=0:MaxPop_F
    for n=0:MaxPop_M
     delta(n+1+(MaxPop_M+1)*m,:)=mean(counts_M'==n & counts_F'==m ); %%basically computing the joint histogram of counts but flattened
    end
end
hist=delta';


N_Fpre=((1:(MaxPop_F+1))-1)'; %vector with possible female occupation numbers in the system
N_Mpre=((1:(MaxPop_M+1))-1)'; %vector with possible male occupation numbers in the system
N_F=reshape(ones(MaxPop_M+1,1)*N_Fpre',MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
N_M=reshape(N_Mpre*ones(1,MaxPop_F+1),MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
NexpAv_F=hist*N_F; %experimental average number of female flies at each bin
NexpAV_M=hist*N_M; %experimental average number of male flies at each bin
NFNMfac=gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)');

%%first iteration is steepest descent with a variable step size
			
presgrad=-TC_logligrad(rootparams , MaxPop_F, MaxPop_M,Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M); %%steepest descent
alpha=TC_linemin(rootparams, presgrad, alpharoot, MaxPop_F, MaxPop_M,Nbins,Tframes,hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M,-presgrad); %%finding an adecuate step size
params=rootparams+alpha*presgrad; %%steepest descent update
params(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
params(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
params(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
%%

pastconjugdir=presgrad; %%updating the previous conjugate direction
pastconjugdir(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
pastconjugdir(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
pastconjugdir(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
pastparams=rootparams; %%updating the previous argument of the loglikelihood
grad=TC_logligrad(params , MaxPop_F, MaxPop_M,Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);

%%now we use nonlinear conjugate gradients to fin the MLE
tolerance=1E-5 %%when the gradient's magnitude is smaller than this, the search stops
counter=0;
tic
%prelike=[TC_logli(params , MaxPop_F, MaxPop_M, Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M)];
prePar=[params'];
			
while sqrt((params-pastparams)'*(params-pastparams))/(MaxPop+2*Nbins)>tolerance %algorith stops when parameters change on average to the 6th decimal
   conjugdir=TC_getconjugdir(params,pastparams,pastconjugdir, MaxPop_F,MaxPop_M,Nbins,Tframes,hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M,grad);  %%preconditioner used is PR, also it can get automatic resets to steepest descent
   conjugdir(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   conjugdir(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   conjugdir(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   alpha=TC_linemin(params,conjugdir,alpharoot, MaxPop_F,MaxPop_M,Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M,grad); %%finding an adecuate step size
   pastparams=params; %%updating the previous argument of the loglikelihood
   pastparams(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   pastparams(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   pastparams(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   params=params+alpha*conjugdir; %%conjugate gradients update
   params(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   params(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   params(2+MaxPop_M)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   pastconjugdir=conjugdir;  %%updating the previous conjugate direction
   %sqrt(TC_logligrad(params ,  MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M)'*TC_logligrad(params , MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M)) %%checking the norm of the gradient
	grad=TC_logligrad(params , MaxPop_F, MaxPop_M,Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
	counter=counter+1;
		%sqrt(grad'*grad);

	%prelike=[prelike;TC_logli(params , MaxPop_F, MaxPop_M, Nbins,Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M)];
	prePar=[prePar;params'];
end
%TC_MLE_like=prelike;
TC_MLE_like=prePar;
counter
toc
if gauge==0
    V_F=params(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
    V_M=params(MaxPop+Nbins+1:end); %%extracting vexation for males
    fr=params(1:MaxPop); %%extracting frustration
    f_TC=reshape(fr,[MaxPop_M+1,MaxPop_F+1]); %reshaped frustration surface
    [stderrors,CovMat_TC]=TC_getCovMat(fr,V_F,V_M, MaxPop_F,MaxPop_M,Nbins,Tframes,0); %getting the covariance matrix for the parameters in the model without gauge transformation 
    ferror_TC=[0;0;stderrors(1:1+MaxPop_M);0;stderrors(2+MaxPop_M:end)]; %%unpacking errors
    V_Ferror=stderrors(MaxPop-2:MaxPop+Nbins-3); %%unpacking errors for V_F
    V_Merror=stderrors(MaxPop+Nbins-2:end);  %%unpacking errors for V_M
%counter
else
   %%index settup for flattened (row major) joint probability distribution
    N_Fpre=((1:(MaxPop_F+1))-1)'; %vector with possible female occupation numbers in the system
    N_Mpre=((1:(MaxPop_M+1))-1)'; %vector with possible male occupation numbers in the system
    N_F=reshape(ones(MaxPop_M+1,1)*N_Fpre',MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
    N_M=reshape(N_Mpre*ones(1,MaxPop_F+1),MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution

    
    V_F=params(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
    V_M=params(MaxPop+Nbins+1:end); %%extracting vexation for males
    fr=params(1:MaxPop)+sum(params(MaxPop+Nbins-1:end))*N_M/Nbins+sum(params(MaxPop+1:MaxPop+Nbins))*N_F/Nbins; %%extracting frustration
    f_TC=reshape(fr,[MaxPop_M+1,MaxPop_F+1]); %%reshaping the frusration to a (MaxPop_F+1)(MaxPop_M+1) matrix 
    [stderrors,CovMat_TC]=TC_getCovMat(fr, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes, 1);  %getting the covariance matrix for the parameters in the model with gauge transformation that sets the average potential to zero
    ferror_TC=reshape([0;stderrors(1:MaxPop-1)],[MaxPop_M+1,MaxPop_F+1]); %%unpacking errors
    V_Ferror=stderrors(MaxPop-1:MaxPop+Nbins-2); %%unpacking errors for V_F
    V_Merror=stderrors(MaxPop+Nbins-1:end);  %%unpacking errors for V_M
end

