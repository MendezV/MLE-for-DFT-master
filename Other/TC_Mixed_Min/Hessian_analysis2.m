

loaddata2;
alpharoot =  0.0010000;
rootparams=rand([MaxPop+2*Nbins,1]);
[f_TC, V_F, V_M, CovMat_TC, ferror_TC, V_Ferror, V_Merror, TC_MLE_like]=TC_MLE(rootparams,alpharoot,counts_M,counts_F,gauge,tau);
plot(TC_MLE_like)
plot(TC_MLE_like(1:end,1:70))
[D,C]=eig(CovMat_TC);
plot(log(abs(sort(diag(C)))),'x')
[xx,yy]=meshgrid((1:MaxPop_F+1)-1,(1:MaxPop_M+1)-1);
scatter3(xx,yy,f_TC, [], f_TC(:),'x')
surf(reshape(D(1:MaxPop,2),[MaxPop_F+1, MaxPop_M+1]))
scatter3(xx,yy,reshape(D(1:MaxPop,1),[MaxPop_F+1, MaxPop_M+1]), [], reshape(D(1:MaxPop,1),[MaxPop_F+1, MaxPop_M+1])(:),'x')
