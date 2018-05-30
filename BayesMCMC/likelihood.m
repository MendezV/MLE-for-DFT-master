function [ Likely ] = likelihood( hist, histerr,f, V, Tframes, MaxPop, Nbins)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%Pobability for the DFT model
z=sum(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
probmat=(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1

Likely=sum(Tframes*sum(((hist-probmat).^2)./histerr,2))/Nbins;

end

