function params=VdW


xdata=linspace(1,10,10);
ydata=4*exp(-3*xdata);
fun = @(x,xdata)x(1)*exp(x(2)*xdata);
x0=[4,-5];
[beta,resnorm,resid,exitflag,x,lambda,J]= lsqcurvefit(fun,x0,xdata,ydata);
ci = nlparci(beta,resid,'jacobian',J)
params=x;
plot(xdata,fun(x,xdata))

end
