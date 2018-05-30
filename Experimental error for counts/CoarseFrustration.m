function [ f_new , ferror_new ] = CoarseFrustration( f_old ,ferror_old, A_old, A_new, fit)

%IN
%-A_old: area of the bins where we extracted the original frustration
%-f_old: MaxPop+1 column vector that corresponds to the frustration for bin
%area A_old
%-ferror_old: MaxPop+1 column vector that corresponds to the one-sigma error on the frustration for bin
%area A_old
%-A_new: area of the bins where we want to determine the frustration
%-fit: type of procedue that will be implemented 0=interpolation,
%1=quadratic fit, 2=VdW fit.


%%performs one of the three possible procedures specified by fit. for
%%fit==0 interpolates the values of the scaled extracted f_old and
%%ferror_old and then determines the value of frustration and its error for
%%integer values at the values of observed density for the new bin area
%%(f_new, ferror_new). Similarly fit==1 and fit==2 do fits to determine the
%%values of frustration according to the fit with error given by the upper
%%and lower bounds in the 95% confidence interval of the fit.


%%OUT
%%-f_new: MaxPop+1 column vector that corresponds to the frustration for bin
%area A_new, goes up to the highest occupation given the large density
%observed in the f_old dataset but for bin area A_new
%%-ferror_new: error on the values of f_new.

MaxPop=size(f_old,1)-1; %%maximum number of flies observed in the dataset that the input frustration is coming from
w=1./(ferror_old+0.0001); %%weights that will be used in the fits correspond to the inverse of the variances on each point of the frustration
x=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system before coarse graining
xx = linspace(0,MaxPop,1000)';


%%rescaling the x axis that corresponds to the number of flies in the new
%%bin size for constant density
x=x*(A_new/A_old); %just the data points that we have
xx=xx*(A_new/A_old);  %% the values that are going to be interpolated
 
if fit==0  %%if we choose the interpolation option
    
    
    ypred = interp1(x,f_old*(A_new/A_old),xx,'spline'); %performing spline interpolation to the rescaled data
    ferr1= interp1(x,f_old*(A_new/A_old)+ferror_old*(A_new/A_old),xx,'spline');  %performing spline interpolation to the rescaled one sigma upper bound
    ferr2= interp1(x,f_old*(A_new/A_old)-ferror_old*(A_new/A_old),xx,'spline');  %performing spline interpolation to the rescaled one sigma lower bound
    ypredci=[ferr2,ferr1];  %formatting
    
  
   

elseif fit==1 %%if we choose the quadratic fit option
    
    start=[0.1,0.2]; %%root for the fit
    modelFun = @(b,x) b(1).*(x.^2)+b(2).*x;  %model that we are going to use to fit the frustration, in this case quadratic
    wnlm = fitnlm(x,f_old*(A_new/A_old),modelFun,start,'Weight',w); %%performing the weighted non-linear fit
    line(xx,predict(wnlm,xx),'color','b') 
    [ypred,ypredci] = predict(wnlm,xx,'Simultaneous',true);   %%fit and confidence bounds, by default gives the 95% confidence interval
   
 
    
elseif fit==2   %%if we choose the VdW fit option
    start=[-0.1,0.2]; %%root for the fit
    modelFun = @(b,x) -x.*log(1-b(1).*x)+b(2).*x; %model that we are going to use to fit the frustration, in this case VdW with a linear term but without constant shift
    wnlm = fitnlm(x,f_old*(A_new/A_old),modelFun,start,'Weight',w);  %%performing the weighted non-linear fit to the appropiatey scaled data
    line(xx,predict(wnlm,xx),'color','b')
    [ypred,ypredci] = predict(wnlm,xx,'Simultaneous',true); %%fit and confidence bounds by default gives the 95% confidence interval
  
elseif fit==3   %%if we choose the VdW fit option
    chi=1.1/2.1;
    [K,E] = ellipke(chi.^2);
    FF=2*E/(pi*sqrt(1-chi.^2));
    start=[-0.1,0.2]; %%root for the fit
    modelFun = @(b,x) (-x.*log(1-b(1).*x)+b(2).*x)*FF; %model that we are going to use to fit the frustration, in this case VdW with a linear term but without constant shift
    wnlm = fitnlm(x,f_old*(A_new/A_old),modelFun,start,'Weight',w);  %%performing the weighted non-linear fit to the appropiatey scaled data
    line(xx,predict(wnlm,xx),'color','b')
    [ypred,ypredci] = predict(wnlm,xx,'Simultaneous',true); %%fit and confidence bounds by default gives the 95% confidence interval
  
end

%%we now extract the frustration at integer values for fitted curve 
MaxPop_new=floor(max(x)); %% maximum number of agents in a bin of size A_new for observed densities in bins of size A_old
N=((1:(MaxPop_new+1))-1)';  %%possible occupations in the bins of size A_new
f_new=zeros(MaxPop_new+1,1); %%declaring the new frustration
ferror_new=zeros(MaxPop_new+1,1);  %%declaring the new errorbars

for i=0:MaxPop_new
    [minimum index] = min(abs(xx-i)); %%value in linspace closer to each of the possible integer occupations
    f_new(i+1)=ypred(index); %%scaled frustration at integer occupation i 
    ferror_new(i+1)=ypredci(index,2)-ypred(index); %%error for the scaled frustration at integer occupation i 
end

%%plotting the result of the above analysis.
plot(x,f_old*(A_new/A_old),'ko', xx,ypred,'b-', xx,ypredci,'r-');
xlabel('x'); ylabel('y');
legend({'original f', 'Weighted fit', '95% Prediction Limits'},'location','SouthEast');
hold on
errorbar(N,f_new,ferror_new)

end

