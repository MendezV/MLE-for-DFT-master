load Mixed.mat;
counts_F=occ_f';
counts_M=occ_m';
counts=counts_M;
Corr;


MaxPop_F=max(max(counts_F)); %maximum observed packing in the system for Females
MaxPop_M=max(max(counts_M)); %maximum observed packing in the system for Males
MaxPop=(MaxPop_M+1)*(MaxPop_F+1); %%including the zeroth fly for both genders
gauge=0;

Nbins=size(counts_M,1); %total number of bins
Tframes=size(counts_M,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)

