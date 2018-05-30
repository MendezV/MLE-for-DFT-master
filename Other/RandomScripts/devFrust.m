function [Nfrust,devfrust,errordev]=devFrust(f,CovMat)
Nmax=size(f,1);
Nfrust=1:(Nmax-1);
devfrust=diff(f);
errordevpre=zeros(size(devfrust));
errordevpre(1)=sqrt(CovMat(1,1));
for i=1:(Nmax-2)
    errordevpre(i+1)=sqrt(CovMat(i,i)+CovMat(i+1,i+1)-2*CovMat(i+1,i));
end
errordev=errordevpre;
end