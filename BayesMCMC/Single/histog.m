function [ histo, histerr, CovMatHist ] = histog( counts ,tau )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
histo=delta';

CovMatHist=zeros(Nbins,MaxPop+1,MaxPop+1);
histerr=sqrt(histo.*(1-histo)./(Tframes+1));

for i=1:Nbins
    CovMatHist(i,:,:)=histo(i,:)'*histo(i,:)./(Tframes+1)-diag(diag(histo(i,:)'*histo(i,:))./(Tframes+1))+diag(histerr(i,:));
end


dlmwrite('histog1.dat', [Nbins,MaxPop+1], 'delimiter',' ');

dlmwrite('histog2.dat', histo, 'delimiter',' ');

system("cat histog1.dat histog2.dat >histog.dat");
system("rm histog1.dat histog2.dat ");

end

