function [ KLDFT, KLNaive, Entropy, JSdiv, JSdivN] = KLDivergence( counts, f ,V ,VNaive )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
histo=delta';


N=((1:(MaxPop+1))-1)'; %%possible bin occupations for the dataset that we are trying to predict
z=sum(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
probmat=(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1


zNaive=sum(exp(-VNaive*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
probmatNaive=(exp(-VNaive*N')./gamma(ones(Nbins,1)*(N+1)'))./(zNaive*ones(1,MaxPop+1)); %size NbinsxMaxPop+1

Uniform=ones(size(probmatNaive))/MaxPop;

KLDFT=zeros(Nbins,1);
KLNaive=zeros(Nbins,1);
Entropy=zeros(Nbins,1);
JSdiv=zeros(Nbins,1);
JSdivN=zeros(Nbins,1);
for i=1:Nbins
    KLpreDFT=histo(i,:).*log(histo(i,:))-histo(i,:).*log(probmat(i,:));
    KLpreDFT(isnan(KLpreDFT))=0;
    KLDFT(i)=sum(KLpreDFT)/log(2);
    
    KLpreNaive=histo(i,:).*log(histo(i,:))-histo(i,:).*log(probmatNaive(i,:));
    KLpreNaive(isnan(KLpreNaive))=0;
    KLNaive(i)=sum(KLpreNaive)/log(2);
    
    preEntro=histo(i,:).*log(histo(i,:));
    preEntro(isnan(preEntro))=0;
    Entropy(i)=-sum(preEntro)/log(2);
    
    a=histo(i,:).*log(histo(i,:))-histo(i,:).*log(0.5*(probmat(i,:)+histo(i,:)));
    a(isnan(a))=0;
    b=probmat(i,:).*log(probmat(i,:))-probmat(i,:).*log(0.5*(probmat(i,:)+histo(i,:)));
    b(isnan(b))=0;
    JSdiv(i)=sum(0.5*(a+b))/log(2);
    
    a=histo(i,:).*log(histo(i,:))-histo(i,:).*log(0.5*(probmatNaive(i,:)+histo(i,:)));
    a(isnan(a))=0;
    b=probmatNaive(i,:).*log(probmatNaive(i,:))-probmatNaive(i,:).*log(0.5*(probmatNaive(i,:)+histo(i,:)));
    b(isnan(b))=0;
    JSdivN(i)=sum(0.5*(a+b));
end

%%entropy of uniform is 165.3750 (maximum entropy), of Naive is 125.6746, of model is 102.1348 and of data hist is  100.7787 KL div
%%of uniform dist is 46.7, of naive is 10.5743 of model is 0.9099
end

