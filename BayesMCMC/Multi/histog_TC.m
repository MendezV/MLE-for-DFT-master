function [ histo ] = histog_TC(counts_M,counts_F,tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



MaxPop_F=max(max(counts_F)); %maximum observed packing in the system for Females
MaxPop_M=max(max(counts_M)); %maximum observed packing in the system for Males
MaxPop=(MaxPop_M+1)*(MaxPop_F+1); %%including the zeroth fly for both genders

Nbins=size(counts_M,1); %total number of bins
Tframes=size(counts_M,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop,Nbins); %%allocating the matrix of kronecker delta functions

for m=0:MaxPop_F
for n=0:MaxPop_M
delta(n+1+(MaxPop_M+1)*m,:)=mean(counts_M'==n & counts_F'==m ); %%basically computing the joint histogram of counts but flattened
end
end
histo=delta';



dlmwrite('histog1.dat', [Nbins,MaxPop], 'delimiter',' ');
dlmwrite('histog2.dat', histo, 'delimiter',' ');

system("cat histog1.dat histog2.dat >histog.dat");
system("rm histog1.dat histog2.dat ");

%%m=int2str(MaxPop_M)
%% c=int2str(MaxPop_F)
%%strcat(c,m)


end

