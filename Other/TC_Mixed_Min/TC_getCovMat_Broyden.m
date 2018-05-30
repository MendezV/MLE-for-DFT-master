function [Hess,CovMat,det1,det2]=TC_getCovMat_Broyden(f, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes,N_F, N_M, NFNMfac)
%%IN
%%-f: a (MaxPop_F+1)(MaxPop_M+1) matrix (MaxPop_F(M) is the maximum number of female (male) flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)
%%-V_F: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for females
%%-V_M: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for males
%%-MaxPop_F:  maximum observed packing in the system for females
%%-MaxPop_M:  maximum observed packing in the system for males
%%-Nbins: total number of bins
%%-Tframes: number of frames
%%-gauge: numerical value, if equal to zero a gauge transformation was
%%performed to the parameters and we have to transform the covariance
%%matrix accordingly


%%Calculates the asymptotic Covariance Matrix as the inverse of the fisher
%%information matrix for our log-likelihood function Same size hessian inverse 


%%OUT
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


%parameters that will be used in the calculation
MaxPop=(MaxPop_F+1)*(MaxPop_M+1); %%including the zeroth fly for both genders


%%frustration sector of the covariance matrix
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2); %%vector that contains the partition function of each bin, size Nbinsx1
probmat=(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac)./(z*ones(1,MaxPop)); %size NbinsxMaxPop+1
			 


				
%%vexation sector of the covariance matrix, ensemble averages
NensAv_F=sum((ones(Nbins,1)*N_F').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2)./z; %%ensemble average of the number of females in our model
NensAv_M=sum((ones(Nbins,1)*N_M').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2)./z; %%ensemble average of the number of males in our model
NsqensAv_F=sum((ones(Nbins,1)*N_F').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2)./z; %%ensemble average of the number of females in our model
NsqensAv_M=sum((ones(Nbins,1)*N_M').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2)./z; %%ensemble average of the number of males in our model
NensAv_M_F=sum(((ones(Nbins,1)*N_M').*(ones(Nbins,1)*N_F')).*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./NFNMfac,2)./z; %%ensemble average of the number of females in our model
	
																 
																 
%%frustration sector of the covariance matrix
HessFF1=Tframes*(diag(sum(probmat,1)')-probmat'*probmat); %%this hessian is the fisher information matrix
HessFF=HessFF1([3:1+MaxPop_M,3+MaxPop_M:end],[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case

%%vexation sector of the covariance matrix
HessV_FV_F=Tframes*diag(NsqensAv_F-NensAv_F.*NensAv_F);
HessV_MV_M=Tframes*diag(NsqensAv_M-NensAv_M.*NensAv_M);
HessV_MV_F=Tframes*diag(NensAv_M_F-NensAv_M.*NensAv_F);


%%mixed sector of V's and f's of the covariance matrix
HessV_FF1=Tframes*probmat.*(ones(Nbins,1)*N_F'- NensAv_F*ones(1,MaxPop));
HessV_MF1=Tframes*probmat.*(ones(Nbins,1)*N_M'- NensAv_M*ones(1,MaxPop));
HessV_FF=HessV_FF1(:,[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case
HessV_MF=HessV_MF1(:,[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case

%%final result

    Hessp=[HessFF,HessV_FF',HessV_MF';HessV_FF,HessV_FV_F,HessV_MV_F;HessV_MF,HessV_MV_F,HessV_MV_M]; %%hessian, equal to the fisher information matrix
	CovMatp=inv(Hessp); %%pseudoinverse, apparently its singular due to gauge invariance
	det1=det(Hessp);
	det2=det(CovMatp);
	Hesspf=zeros([MaxPop+2*Nbins,MaxPop+2*Nbins]);
    CovMatpf=zeros([MaxPop+2*Nbins,MaxPop+2*Nbins]);
	
	CovMatp1=CovMatp(1:MaxPop_M-1,1:MaxPop_M-1);
	CovMatp2=CovMatp(1:MaxPop_M-1,MaxPop_M:end);
	CovMatp3=CovMatp(MaxPop_M:end,1:MaxPop_M-1);
	CovMatp4=CovMatp(MaxPop_M:end,MaxPop_M:end);
							
							
	Hessp1=Hessp(1:MaxPop_M-1,1:MaxPop_M-1);
	Hessp2=Hessp(1:MaxPop_M-1,MaxPop_M:end);
	Hessp3=Hessp(MaxPop_M:end,1:MaxPop_M-1);
	Hessp4=Hessp(MaxPop_M:end,MaxPop_M:end);
							
	CovMatpf(3:MaxPop_M+1,3:MaxPop_M+1)=CovMatp1;
	CovMatpf(3:MaxPop_M+1,MaxPop_M+3:end)=CovMatp2;
	CovMatpf(MaxPop_M+3:end,3:MaxPop_M+1)=CovMatp3;
	CovMatpf(MaxPop_M+3:end,MaxPop_M+3:end)=CovMatp4;
							
	Hesspf(3:MaxPop_M+1,3:MaxPop_M+1)=Hessp1;
	Hesspf(3:MaxPop_M+1,MaxPop_M+3:end)=Hessp2;
	Hesspf(MaxPop_M+3:end,3:MaxPop_M+1)=Hessp3;
    Hesspf(MaxPop_M+3:end,MaxPop_M+3:end)=Hessp4;
							
	Hess=Hesspf;
	CovMat=CovMatpf;



end
