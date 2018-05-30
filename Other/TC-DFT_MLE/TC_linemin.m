function alpha=TC_linemin(params,conjugdir,alpharoot, MaxPop_F, MaxPop_M, Nbins, Tframes, hist)

%%IN
%%-params: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins, corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-conjugdir: a vector of size (MaxPop_F+1)(MaxPop_M+1) + 2*Nbins  corresponding to the
%%new search direction in the algorithm
%%-alpharoot: scalar that corresponds to the root that will be used in
%%the linesearch algorithm
%%-MaxPop_F:  maximum observed packing in the system for females
%%-MaxPop_M:  maximum observed packing in the system for males
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%performs a backtracking search to determine a sufficiently good step size
%%to guaratee that we are indeed minimizing the function. Sufficiently good
%%is determined by the Armijo-Goldstein condition

%%OUT
%%-alpha: optimal step size along the search direction according to the Armijo-Goldstein condition

tau=0.7; %%control parameter should be between 0 and 1
c1=0.0001; %%control parameter should be between 0 and 1
c2=0.1;  %%control parameter should be between 0 and 1
grad=TC_logligrad(params ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist); %%gradient of the function to be minimized
m=conjugdir'*grad/sqrt(conjugdir'*conjugdir); %%local slope of the function of alpha  along the search direction
t=-c1*m; %parameter that serves as lower bound in the Armijo?Goldstein condition 
alpha=alpharoot; %%initializing alpha for the linesearch

diff=TC_logli(params , MaxPop_F, MaxPop_M, Nbins,Tframes, hist)-TC_logli(params+alpha*conjugdir , MaxPop_F, MaxPop_M, Nbins,Tframes, hist); %%difference that will be minimized 
grad2=TC_logligrad(params+ alpha*conjugdir, MaxPop_F, MaxPop_M, Nbins, Tframes, hist);
%%updating alpha until the Armijo-Goldstein and the sstrong Wolfe conditions are fulfilled, then we
%%are satisfied with the value of alpha that will multiply the conjugate
%%direction in the algorithm

while diff<alpha*t && abs(conjugdir'*grad2)>c2*abs(conjugdir'*grad)
    alpha=alpha*tau;
    diff=TC_logli(params , MaxPop_F, MaxPop_M, Nbins,Tframes, hist)-TC_logli(params+alpha*conjugdir , MaxPop_F, MaxPop_M, Nbins,Tframes, hist); %%difference that will be minimized 
    grad2=TC_logligrad(params+ alpha*conjugdir, MaxPop_F, MaxPop_M, Nbins, Tframes, hist);
    
end

end

