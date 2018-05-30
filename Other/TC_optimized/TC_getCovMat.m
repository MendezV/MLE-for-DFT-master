function [stderrors,CovMat]=TC_getCovMat(f, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes, gauge)

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
%%information matrix for our log-likelihood function


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

%%index settup for flattened (row major) joint probability distribution
N_Fpre=((1:(MaxPop_F+1))-1)'; %vector with possible female occupation numbers in the system
N_Mpre=((1:(MaxPop_M+1))-1)'; %vector with possible male occupation numbers in the system
N_F=reshape(ones(MaxPop_M+1,1)*N_Fpre',MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
N_M=reshape(N_Mpre*ones(1,MaxPop_F+1),MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution


%%frustration sector of the covariance matrix
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
probmat=(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')))./(z*ones(1,MaxPop)); %size NbinsxMaxPop+1
HessFF=Tframes*(diag(sum(probmat,1)')-probmat'*probmat); %%this hessian is the fisher information matrix 
HessFF=HessFF([3:1+MaxPop_M,3+MaxPop_M:end],[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case

%%vexation sector of the covariance matrix
NensAv_F=sum((ones(Nbins,1)*N_F').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NsqensAv_F=sum((ones(Nbins,1)*(N_F.^2)').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
HessV_FV_F=Tframes*diag(NsqensAv_F-NensAv_F.*NensAv_F);

NensAv_M=sum((ones(Nbins,1)*N_M').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NsqensAv_M=sum((ones(Nbins,1)*(N_M.^2)').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
HessV_MV_M=Tframes*diag(NsqensAv_M-NensAv_M.*NensAv_M);

NensAv_M_F=sum(((ones(Nbins,1)*N_M').*(ones(Nbins,1)*N_F')).*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*f')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
HessV_MV_F=Tframes*diag(NensAv_M_F-NensAv_M.*NensAv_F);


%%mixed sector of V's and f's of the covariance matrix
HessV_FF=Tframes*probmat.*(ones(Nbins,1)*N_F'- NensAv_F*ones(1,MaxPop));
HessV_FF=HessV_FF(:,[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case
HessV_MF=Tframes*probmat.*(ones(Nbins,1)*N_M'- NensAv_M*ones(1,MaxPop));
HessV_MF=HessV_MF(:,[3:1+MaxPop_M,3+MaxPop_M:end]);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case


%%final result
if gauge==0
    Hess=[HessFF,HessV_FF',HessV_MF';HessV_FF,HessV_FV_F,HessV_MV_F;HessV_MF,HessV_MV_F,HessV_MV_M]; %%hessian, equal to the fisher information matrix
    %CovMat=inv(Hess); %%pseudoinverse, apparently its singular due to gauge invariance
	CovMat=pinv(Hess);
    stderrors=sqrt(diag(CovMat)); %%the factor of two comes in the taylor expansion for the gaussian approximation

else
    Hess=[HessFF,HessV_FF',HessV_MF';HessV_FF,HessV_FV_F,HessV_MV_F;HessV_MF,HessV_MV_F,HessV_MV_M]; %%hessian, equal to the fisher information matrix
    CovMatp=inv(Hess); %%pseudoinverse, apparently its singular due to gauge invariance
    gaugeTMat=eye(size(Hess))+[zeros(size(HessFF)),N_F([3:1+MaxPop_M,3+MaxPop_M:end])*ones(1,Nbins)./Nbins,(N_M([3:1+MaxPop_M,3+MaxPop_M:end]))*ones(1,Nbins)./Nbins;zeros(size(HessV_FF)),-ones(size(HessV_FV_F))./Nbins,zeros(size(HessV_MV_F));zeros(size(HessV_MF)),zeros(size(HessV_MV_F)),-ones(size(HessV_MV_M))./Nbins];
    CovMatg=gaugeTMat*CovMatp*gaugeTMat'; %%the gauge transformation to set the average potential is just a linear transformation on the vector of parameters, this is the induced transformation on the covariance matrix
    
    %sectors of the non-gauge transformed covariance matrix
    CovMatFV_F=CovMatp(1:MaxPop-3,MaxPop-2:MaxPop-2+Nbins-1);
    CovMatFV_M=CovMatp(1:MaxPop-3,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1);
    CovMatV_FV_F=CovMatp(MaxPop-2:MaxPop-2+Nbins-1,MaxPop-2:MaxPop-2+Nbins-1);
    CovMatV_MV_M=CovMatp(MaxPop-2+Nbins:MaxPop-2+2*Nbins-1,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1);
    CovMatV_MV_F=CovMatp(MaxPop-2:MaxPop-2+Nbins-1,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1);
    %CovMatFF=CovMatp(1:MaxPop-3,1:MaxPop-3);
    
    
    %%the new elements that are going to be added to the covariance matrix
    %%and are calculated to first order for standard error propagation
    %%taking into account correlations between the gauge transformed
    %%parameters. 10==1 male 0 females 01== 0 males 1 female.
    varf10=sum(sum(CovMatV_MV_M))./(Nbins^2); 
    varf01=sum(sum(CovMatV_FV_F))./(Nbins^2);
    covf01f10=sum(sum(CovMatV_MV_F))./(Nbins.^2);
    covf01V_F=sum(CovMatV_FV_F)./Nbins+sum(sum(CovMatV_FV_F))./(Nbins^2);
    covf01V_M=sum(CovMatV_MV_F)./Nbins+sum(sum(CovMatV_MV_F))./(Nbins^2);
    covf01F=sum(CovMatFV_F,2)./Nbins+N_M([3:1+MaxPop_M,3+MaxPop_M:end])*varf10/(Nbins^2)+N_F([3:1+MaxPop_M,3+MaxPop_M:end])*varf01/(Nbins^2);
    covf10V_F=sum(CovMatV_MV_F)./Nbins+sum(sum(CovMatV_MV_F))./(Nbins^2);
    covf10V_M=sum(CovMatV_MV_M)./Nbins+sum(sum(CovMatV_MV_M))./(Nbins^2);
    covf10F=sum(CovMatFV_M,2)./Nbins+N_M([3:1+MaxPop_M,3+MaxPop_M:end])*varf10/(Nbins^2)+N_F([3:1+MaxPop_M,3+MaxPop_M:end])*varf01/(Nbins^2);
    
    %%partition of the matrix to include the new sectors with covariances
    %%and variances of the redundant parameters 
    CovMatg1=CovMatg(1:MaxPop_M-1,1:MaxPop_M-1);
    CovMatg2=CovMatg(1:MaxPop_M-1,MaxPop_M:end);
    CovMatg3=CovMatg(MaxPop_M:end,1:MaxPop_M-1);
    CovMatg4=CovMatg(MaxPop_M:end,MaxPop_M:end);
    
    %%transpossing what needs to be transposed
    covf01V_F=covf01V_F';
    covf01V_M=covf01V_M';
    covf10V_F=covf10V_F';
    covf10V_M=covf10V_M';
   
    
    %%assembling the new covariance matrix
    CovMat=[varf01,covf01F(1:MaxPop_M-1)',covf01f10,covf01F(MaxPop_M:end)', covf01V_F', covf01V_M';covf01F(1:MaxPop_M-1),CovMatg1,covf10F(1:MaxPop_M-1),CovMatg2;covf01f10,covf10F(1:MaxPop_M-1)',varf10,covf10F(MaxPop_M:end)',covf10V_F', covf10V_M';[covf01F(MaxPop_M:end); covf01V_F; covf01V_M],CovMatg3,[covf10F(MaxPop_M:end); covf10V_F; covf10V_M],CovMatg4]; %% covariance matrix with the f(1) uncertainties and covariances added to it, f(1) variance includes a sum of al the variances and coariances of the v's divided by the number of bins squared, the other new sectors are the covariances with the other parameters, which can be easily derived by considering the linearity of the average and the definition of the covariance
    stderrors=sqrt(diag(CovMat)); %diagonal of the matrix above corresponds to the variances of the parameters

end


end
