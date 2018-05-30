function [predictedAv_M, predictedAv_F, errorbars_F, errorbars_M]=TC_predav(f, V_F, V_M, CovMatF, CovMatV_F, CovMatV_M, Nflies_M, Nflies_F,sameexperiment,gauge)

%%IN:
%%-f: a (MaxPop_F+1)(MaxPop_M+1) matrix (MaxPop_F(M) is the maximum number of female (male) flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)
%%-V_F: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for females
%%-V_M: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average
%%potential is set to zero) for males
%%-CovMatV_M a ((MaxPop_F+1)(MaxPop_M+1)-3+2*Nbins) square covariance matrix obtained from
%%getCovMatF: it has less dimensions than expected from the number of parameters in the model due to the lack of error
%%in the gauge fixed parameters. Can also be a 2*Nbins square matrix.
%%-CovMatV_F a ((MaxPop_F+1)(MaxPop_M+1)-3+2*Nbins) square covariance matrix obtained from
%%getCovMatF: it has less dimensions than expected from the number of parameters in the model due to the lack of error
%%in the gauge fixed parameters. Can also be a 2*Nbins square matrix.
%%-CovMatF a ((MaxPop_F+1)(MaxPop_M+1)-3+2*Nbins) square covariance matrix obtained from
%%getCovMatF: it has less dimensions than expected from the number of parameters in the model due to the lack of error
%%in the gauge fixed parameters.
%%-Nflies_F: Number of female flies in the sistem that we are trying to predict
%%-Nflies_M: Number of male flies in the sistem that we are trying to predict
%%-sameexperiment: int, if it is 1 we take into account correlations between parameters V and F otherwise they will not be correlated 
%%-gauge numerical value, if 0 no gauge transformation is applied to the
%%parameters and f(0)=f(1)=0. else the gauge is set so that the average
%%potential is zero and f(1) corresponds to the sum of the previous values
%%for the potential


%%predicts the average number of flies in each bin using the DFT model and
%%then determines the uncertainties taking into account correlations
%%between the parameters in the model

%%OUT:
%%-predictedAv_F a vector of size Nbins (the number of bins in the system) in
%%which each entry is the average number of female flies as predicted from the DFT
%%model implemented with a binary search that determines the chemical
%%potential mu_F that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)
%%-errorbars_F a vector of size Nbins (the number of bins in the system) in
%%which each entry is the one sigma error on the average number of female flies as predicted from the DFT
%%model for each bin(vector is sorted in the same order as the occupations matrix first dimension)
%%-predictedAv_M a vector of size Nbins (the number of bins in the system) in
%%which each entry is the average number of male flies as predicted from the DFT
%%model implemented with a binary search that determines the chemical
%%potential mu_M that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)
%%-errorbars_M a vector of size Nbins (the number of bins in the system) in
%%which each entry is the one sigma error on the average number of male flies as predicted from the DFT
%%model for each bin(vector is sorted in the same order as the occupations matrix first dimension)



MaxPop_M=size(f,1)-1; % %maximum observed packing for males in the system where F was extracted (must be greater or equal to the maximum observed packing in the data we are about to predict)
MaxPop_F=size(f,2)-1; %maximum observed packing for females in the system where F was extracted (must be greater or equal to the maximum observed packing in the data we are about to predict)
Nbins=size(V_F,1); %total number of bins in the data we are about to predict
MaxPop=(MaxPop_M+1)*(MaxPop_F+1); %%including the zeroth fly for both genders

%%index settup for flattened (row major) joint probability distribution
N_Fpre=((1:(MaxPop_F+1))-1)'; %vector with possible female occupation numbers in the system
N_Mpre=((1:(MaxPop_M+1))-1)'; %vector with possible male occupation numbers in the system
N_F=reshape(ones(MaxPop_M+1,1)*N_Fpre',MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution
N_M=reshape(N_Mpre*ones(1,MaxPop_F+1),MaxPop,1); %%reshaping so that dimensions match with the flattened joint probability distribution

%% first we make the prediction for females and males

NmaxSteps=200;    %number of steps before giving up on the search
tolerance= 1e-7;  %tolerance in the difference between the actual number of flies and the one that comes up with the guessed mu
counter=0;
mu_F=max(V_F); %%convenient root for the search algorithm
mu_M=max(V_M); %%convenient root for the search algorithm
NfliesGuess_F=0.0;
NfliesGuess_M=0.0;
a_F=mu_F+5;  %guess for the highest lower bound for the chemical potential
b_F=mu_F-5;   %guess for the lowest higher bound for the chemical potential
a_M=mu_M+5;  %guess for the highest lower bound for the chemical potential
b_M=mu_M-5;   %guess for the lowest higher bound for the chemical potential
remu_F=0;
remu_M=0;

F=reshape(f,[MaxPop,1]);
while counter<NmaxSteps
   
    mu_F=(a_F+b_F)/2;
    mu_M=(a_M+b_M)/2;
    z=sum(exp(-(V_F-mu_F)*N_F'-(V_M-mu_M)*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
    predictedAv_F=sum((ones(Nbins,1)*N_F').*exp(-(V_F-mu_F)*N_F'-(V_M-mu_M)*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
    predictedAv_M=sum((ones(Nbins,1)*N_M').*exp(-(V_F-mu_F)*N_F'-(V_M-mu_M)*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
    NfliesGuess_F=sum(predictedAv_F); %% total number of female flies for this iteration of mu_F
    NfliesGuess_M=sum(predictedAv_M); %% total number of male flies for this iteration of mu_F
  
   if abs(NfliesGuess_F-Nflies_F)<tolerance && abs(NfliesGuess_M-Nflies_M)<tolerance%% if the number of flies is right we have the "real mu" and break the while
       remu_F=mu_F; %%real mu in case we need it or something
       remu_M=mu_M; %%real mu in case we need it or something
       break
   end
   
   counter=counter+1; %% if the number of flies is not right counter goes up and we modify accordingly for the binary search
   if(NfliesGuess_F<Nflies_F)
		b_F=mu_F;
		
   else
		a_F=mu_F;
   end
   if(NfliesGuess_M<Nflies_M)
		b_M=mu_M;
		
   else
		a_M=mu_M;
   end
    
end


%% then we calculate the error on that prediction by standard propagation of errors
%%these procedurres were implemented in a similar way to calculate elements
%%of the covariance matrix but we have to repeat them due to the presence
%%of the chemical potential we have just introduced


%frist some useful averages
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
NsqensAv_F=sum((ones(Nbins,1)*(N_F.^2)').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NensAv_F=sum((ones(Nbins,1)*N_F').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NensAv_M=sum((ones(Nbins,1)*N_M').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NsqensAv_M=sum((ones(Nbins,1)*(N_M.^2)').*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model
NensAv_M_F=sum(((ones(Nbins,1)*N_M').*(ones(Nbins,1)*N_F')).*exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2)./z; %%ensemble average of the number of females in our model

%% FOR THE FEMALE AVERAGE%%%%%%%%

%%with respect to V_F
partialV_F=-(NsqensAv_F-NensAv_F.*NensAv_F);
stderrors=sqrt(diag(CovMatV_F)); %%errors taken from the diagonal of the matrix CovMatV_F due to the separability of the distriution the covariances between vexations don't enter in the calculation
V_Fexerror=stderrors(end-2*Nbins+1:end-Nbins);

%%with respect to V_M
partialV_M=-(NensAv_M_F-NensAv_M.*NensAv_F);
stderrors=sqrt(diag(CovMatV_M)); %%errors taken from the diagonal of the matrix CovMatV_M due to the separability of the distriution the covariances between vexations don't enter in the calculation
V_Mexerror=stderrors(end-Nbins+1:end);

%%correlations between values of V do not appear as the partial derivative with respect to
%%vexations in other bins vanishes, the average in each bin doesnt depend on the
%%vexation in other bins

%%with respect to f
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
probmat=(exp(-(V_F-remu_F)*N_F'-(V_M-remu_M)*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')))./(z*ones(1,MaxPop)); %size NbinsxMaxPop+1
partialF=-probmat.*(ones(Nbins,1)*N_F'- NensAv_F*ones(1,MaxPop));

if gauge==0 %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    partialF=partialF(:,[3:1+MaxPop_M,3+MaxPop_M:end]); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case
    
    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters

    partialFpartialF=zeros(Nbins,MaxPop-3,MaxPop-3);
    corrFF=CovMatF(1:MaxPop-3,1:MaxPop-3);
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

else %%condition if the gauge is fixed so that the average potential is zero
    
    partialF=partialF(:,2:end); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case

    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters, but we also have to take into account the
    %%correlations with the value of f(1)

    partialFpartialF=zeros(Nbins,MaxPop-1,MaxPop-1);
    corrFF=CovMatF(1:MaxPop-1,1:MaxPop-1); %%multiplying off diagonal elemnts by 2 right off the bat
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

end

%%the following analysis just applies to same experiment case

if gauge==0   %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    %%for both taking into account correlated errors
    partialFpartialV_F=partialF.*(partialV_F*ones(1,MaxPop-3)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_F=CovMatF(1:MaxPop-3,MaxPop-2:MaxPop-2+Nbins-1)'; %%add two to each index for the non gauge fixed case
    partialFpartialV_M=partialF.*(partialV_M*ones(1,MaxPop-3)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_M=CovMatF(1:MaxPop-3,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1)';
    partialV_MpartialV_F=partialV_M.*partialV_F;
    CovMatV_MV_F=diag(CovMatF(MaxPop-2:MaxPop-2+Nbins-1,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1));


else  %%condition if the gauge is fixed so that the average potential is zero
   
    %%for both taking into account correlated errors and including
    %%correlations with f(0,1) and f(1,0)
    partialFpartialV_F=partialF.*(partialV_F*ones(1,MaxPop-1)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_F=CovMatF(1:MaxPop-1,MaxPop:MaxPop+Nbins-1)'; 
    partialFpartialV_M=partialF.*(partialV_M*ones(1,MaxPop-1)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_M=CovMatF(1:MaxPop-1,MaxPop+Nbins:MaxPop+2*Nbins-1)';
    partialV_MpartialV_F=partialV_M.*partialV_F;
    CovMatV_MV_F=diag(CovMatF(MaxPop:MaxPop+Nbins-1,MaxPop+Nbins:MaxPop+2*Nbins-1));

end



%%now we calculate the error bars with standard propagation of errors to
%%first order taking into account correlations 

if sameexperiment==1 %%condition if both the V's and F come from the same experiment
    errorbars_F=sqrt((partialV_M.^2).*(V_Mexerror.^2)+(partialV_F.^2).*(V_Fexerror.^2)+ferrorsq+sum(partialFpartialV_F.*CovMatFV_F,2)+sum(partialFpartialV_M.*CovMatFV_M,2)+2*partialV_MpartialV_F.*CovMatV_MV_F);
else  %%condition if both the V's and F come from three different experiments
    errorbars_F=sqrt((partialV_F.^2).*(V_Fexerror.^2)+(partialV_M.^2).*(V_Mexerror.^2)+ferrorsq);
end

%%%%%%%%%%



%% FOR THE MALE AVERAGE

%%with respect to V_M
partialV_M=-(NsqensAv_M-NensAv_M.*NensAv_M);
stderrors=sqrt(diag(CovMatV_M)); %%errors taken from the diagonal of the matrix CovMatV_F due to the separability of the distriution the covariances between vexations don't enter in the calculation
V_Mexerror=stderrors(end-Nbins+1:end);

%%with respect to V_F
partialV_F=-(NensAv_M_F-NensAv_M.*NensAv_F);
stderrors=sqrt(diag(CovMatV_M)); %%errors taken from the diagonal of the matrix CovMatV_M due to the separability of the distriution the covariances between vexations don't enter in the calculation
V_Fexerror=stderrors(end-2*Nbins+1:end-Nbins);

%%correlations between values of V do not appear as the partial derivative with respect to
%%vexations in other bins vanishes, the average in each bin doesnt depend on the
%%vexation in other bins

%%with respect to f
z=sum(exp(-V_F*N_F'-V_M*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')),2); %%vector that contains the partition function of each bin, size Nbinsx1
probmat=(exp(-(V_F-remu_F)*N_F'-(V_M-remu_M)*N_M'-ones(Nbins,1)*F')./(gamma(ones(Nbins,1)*(N_F+1)').*gamma(ones(Nbins,1)*(N_M+1)')))./(z*ones(1,MaxPop)); %size NbinsxMaxPop+1
partialF=-probmat.*(ones(Nbins,1)*N_M'- NensAv_M*ones(1,MaxPop));

if gauge==0 %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    partialF=partialF(:,[3:1+MaxPop_M,3+MaxPop_M:end]); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case
    
    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters

    partialFpartialF=zeros(Nbins,MaxPop-3,MaxPop-3);
    corrFF=CovMatF(1:MaxPop-3,1:MaxPop-3);
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

else %%condition if the gauge is fixed so that the average potential is zero
    
    partialF=partialF(:,2:end); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case

    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters, but we also have to take into account the
    %%correlations with the value of f(1)

    partialFpartialF=zeros(Nbins,MaxPop-1,MaxPop-1);
    corrFF=CovMatF(1:MaxPop-1,1:MaxPop-1); %%multiplying off diagonal elemnts by 2 right off the bat
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

end

%%the following analysis just applies to same experiment case

if gauge==0   %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    %%for both taking into account correlated errors
    partialFpartialV_F=partialF.*(partialV_F*ones(1,MaxPop-3)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_F=CovMatF(1:MaxPop-3,MaxPop-2:MaxPop-2+Nbins-1)'; 
    partialFpartialV_M=partialF.*(partialV_M*ones(1,MaxPop-3)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_M=CovMatF(1:MaxPop-3,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1)';
    partialV_MpartialV_F=partialV_M.*partialV_F;
    CovMatV_MV_F=diag(CovMatF(MaxPop-2:MaxPop-2+Nbins-1,MaxPop-2+Nbins:MaxPop-2+2*Nbins-1));


else  %%condition if the gauge is fixed so that the average potential is zero
   
    %%for both taking into account correlated errors and including
    %%correlations with f(0,1) and f(1,0)
    partialFpartialV_F=partialF.*(partialV_F*ones(1,MaxPop-1)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_F=CovMatF(1:MaxPop-1,MaxPop:MaxPop+Nbins-1)'; %%add two to each index for the non gauge fixed case
    partialFpartialV_M=partialF.*(partialV_M*ones(1,MaxPop-1)); %%formula is the multiplication of the partial derivatives times the covariance
    CovMatFV_M=CovMatF(1:MaxPop-1,MaxPop+Nbins:MaxPop+2*Nbins-1)';
    partialV_MpartialV_F=partialV_M.*partialV_F;
    CovMatV_MV_F=diag(CovMatF(MaxPop:MaxPop+Nbins-1,MaxPop+Nbins:MaxPop+2*Nbins-1));

    

end

%%now we calculate the error bars with standard propagation of errors to
%%first order taking into account correlations 

if sameexperiment==1 %%condition if both V's and F come from the same experiment
    errorbars_M=sqrt((partialV_M.^2).*(V_Mexerror.^2)+(partialV_F.^2).*(V_Fexerror.^2)+ferrorsq+sum(partialFpartialV_F.*CovMatFV_F,2)+sum(partialFpartialV_M.*CovMatFV_M,2)+2*partialV_MpartialV_F.*CovMatV_MV_F);
else  %%condition if both V's and F come from different experiments
    errorbars_M=sqrt((partialV_F.^2).*(V_Fexerror.^2)+(partialV_M.^2).*(V_Mexerror.^2)+ferrorsq);
end

%%%%%%%%%%

 

end