function [paramss,Hess, CovMat,TC_MLE_like]=Broyden_TC(params, pastparams , MaxPop_F, MaxPop_M, Nbins, Tframes, hist, prePar, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M)
 
 %%%
%%%%%% everything here fails due to assumptions of Broyden and newton methods!!!!
MaxPop=(MaxPop_F+1)*(MaxPop_M+1); %%including the zeroth fly for both genders

 V_F=params(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
 V_M=params(MaxPop+Nbins+1:end); %%extracting vexation for males
 fr=params(1:MaxPop); %%extracting frustration
 [Jpas,Jinvpas]=TC_getCovMat_Broyden(fr, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes,N_F, N_M, NFNMfac);


 xn=params
 xnpas=pastparams;
 fnpas=TC_logligrad(xnpas ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
%{
 for i=1:20
    fn=TC_logligrad(xn ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);

	V_F=xn(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
	V_M=xn(MaxPop+Nbins+1:end); %%extracting vexation for males
	fr=xn(1:MaxPop); %%extracting frustration
	[J_teo,Jinv_teo]=TC_getCovMat_Broyden(fr, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes,N_F, N_M, NFNMfac);

    deltaxn=xn-xnpas;
    deltafn=fn-fnpas;
 
    deltaxnnorm=sum(deltaxn.^2);
	J=Jpas+(1/deltaxnnorm)*(deltafn-Jpas*deltaxn)*deltaxn';
    normJinv=deltaxn'*(Jinvpas*deltafn);
	Jinv=Jinvpas+((deltaxn-Jinvpas*deltafn)*deltaxn')*Jinvpas./normJinv;
	delta=-Jinv*fn;
				  
	delta(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
	delta(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
	delta(2+MaxPop_M)=0;
	xnfut=xn+delta
	sqrt(sum(fn.^2));
 
 	xnpas=xn;
 	xn=xnfut;
 	fnpas=fn;
 	Jinvpas=Jinv;
 	Jpas=J;
							   
	prePar=[prePar;xn'];
 end

	%}

fn=TC_logligrad(xn ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
sqrt(sum(fn.^2))
V_F=xn(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
V_M=xn(MaxPop+Nbins+1:end); %%extracting vexation for males
fr=xn(1:MaxPop); %%extracting frustration
			
[J,Jinv]=TC_getCovMat_Broyden(fr, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes,N_F, N_M, NFNMfac);
%{
 for i=1:5
			V_F=xn(MaxPop+1:MaxPop+Nbins); %%extracting vexation for females
			V_M=xn(MaxPop+Nbins+1:end); %%extracting vexation for males
			fr=xn(1:MaxPop); %%extracting frustration
			[J,Jinv,det1,det2]=TC_getCovMat_Broyden(fr, V_F, V_M, MaxPop_F, MaxPop_M, Nbins, Tframes,N_F, N_M, NFNMfac);

			fn=TC_logligrad(xn ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
			delta=-Jinv*fn;
			xn=xn+delta
			det1
			det2
			
			sqrt(sum(fn.^2))
			
			
		    prePar=[prePar;xn'];
end
					%}
					
for i=1:1000
					fn=TC_logligrad(xn ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
					fnpas=TC_logligrad(xnpas ,MaxPop_F, MaxPop_M, Nbins, Tframes, hist, N_F, N_M, NFNMfac, NexpAv_F, NexpAV_M);
					deltaxn=xn-xnpas;
					deltafn=fn-fnpas;

					delta=-(deltaxn./(deltafn)).*fn;
					xnfut=xn+delta;
					
					xnpas=xn;
					xn=xnfut;
					
					prePar=[prePar;xn'];
					sqrt(sum(fn.^2))
end

 paramss=xn;
 Hess=J;
 CovMat=Jinv;
 TC_MLE_like=prePar;
			
end
