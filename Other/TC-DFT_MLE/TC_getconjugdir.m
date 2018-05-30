function conjugdir=TC_getconjugdir(params,pastparams,pastconjugdir, MaxPop_F, MaxPop_M, Nbins, Tframes, hist)

%%IN
%%-params: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins, corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastparams corresponds to the position for past
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastconjugdir: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins  corresponding to the
%%last search direction in the algorithm
%%-MaxPop_F:  maximum observed packing in the system for females
%%-MaxPop_M:  maximum observed packing in the system for males
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%Gets the search direction in the non-linear conjugate gradients minimization algorithm by implementing a preconditioner given by the
%%Polak-Ribiere (Fletcher-Reeves) formula with a modification to generate automatic resets to
%%steepest descent if necesary
%%


%%OUT
%%-conjugdir: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins  corresponding to the
%%new search direction in the algorithm

pastgrad=TC_logligrad(pastparams ,  MaxPop_F, MaxPop_M, Nbins, Tframes, hist); %gradient in the last iteration
grad=TC_logligrad(params , MaxPop_F, MaxPop_M, Nbins, Tframes, hist); %gradient in this iteration
%betaPR=grad'*(grad-pastgrad)/(pastgrad'*pastgrad); %Polak-Ribiere 
%beta=max([0,betaPR]); %%automatic reset for steepest decent
betaFR=grad'*grad/(pastgrad'*pastgrad); %Fletcher-Reeves
beta=betaFR;
conjugdir=-grad+beta*pastconjugdir; %%conjugate direction in the algorithm
end