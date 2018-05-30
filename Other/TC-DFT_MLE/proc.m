
n=counts';
fnmin=0;

% ===============================================================
% Autocorrelation function (using FFT), averaged over boxes
c=(mean( abs(ifft(abs(fft(n)).^2))/size(n,1) - ones(size(n,1),1)*mean(n).^2 ,2));
figure(1); plot(c(1:100),'rx')

lit=[1:20]'; % Yunus data seems to show two time scales??
p=polyfit(lit,log(c(lit)),1); 
fprintf("Correlation time (in frames)...\n");
tau=-1/p(1)
fprintf("\n");


%lit=[1:40]';
%lit=[1:40]';
%cdat=[lit-1, polyval(p,lit), log(c(lit))];
%figure(2); plot(cdat(:,1),cdat(:,2),'b-',cdat(:,1),cdat(:,3),'rx')
%save -ascii cdat cdat
%save -ascii c c

% ===============================================================
% Autocorrelation function using dn=n-<n>
dn=n-ones(size(n,1),1)*mean(n); % Compute fluctuations
%Below didn't work as expected for some reason...
% Autocorrelation function (using FFT), averaged over boxes
%c=mean( abs(ifft(abs(fft(dn)).^2)) , 2)/size(dn,1);
%
%lit=[1:30]';
%p=polyfit(lit,log(c(lit)),1); tau=-1/p(1)
%
%lit=[1:40]';
%cdat=[lit-1, polyval(p,lit), log(c(lit))];
%plot(cdat(:,1),cdat(:,2),'b-',cdat(:,1),cdat(:,3),'rx')
%save -ascii cdat cdat
%save -ascii c c

% ===============================================================
% Spatial correlation function
g=dn'*dn;
save -ascii g g

% Make surf plot
%figure(3); surf(g); % In plot window
%gs=smooth(smooth(smooth(smooth(g)))); ppm("g.ppm",gs,gs,gs); system("pnmtopng g.ppm > g.png"); % png gray scale graphic


% ===============================================================
% Computing histograms for each bin
nmx=max(max(n));

h=[];
for cn=0:nmx
  h=[h; sum(n==cn)];
end

nfact=gamma([1:nmx+1]'*ones(1,size(h,2)));
ns=[0:nmx]'*ones(1,size(h,2));

% Basic plots comparing -log(h) -log(n! h) -- showing n! is "correct"
%figure(2); plot(-log(h))
%figure(3); plot(-log(h.*nfact))
%figure(1);

% Probabilies, fractional 1/sigma^2, and -ln (n! p)
p=h/size(n,1);
sis=h/tau;
rho=-log(nfact.*p); rho(sis==0)=0; % Handle Inf issue

% Subsampled matrices for fitting reduced set of f_n's...
sisc=sis(1+fnmin:nmx+1,:);
rhoc=rho(1+fnmin:nmx+1,:);

% Least squares fit (rho, sis, etc., matrices trasposed compared to derivation on board)

% ===============================================================
% Define LHS matrix blocks, and contruct full matrix
ff=diag(sum(sisc,2));
fA=sisc;
fV=diag([fnmin:nmx])*sisc;

Af=sisc';
AA=diag(sum(sis));
AV=diag(sum(ns.*sis));

Vf=sisc'*diag([fnmin:nmx]);
VA=diag(sum(sis.*ns));
VV=diag(sum(ns.^2.*sis));

L=[ff,fA,fV; Af,AA,AV; Vf,VA,VV];

%Define RHS
rhs_f=sum(sisc.*rhoc,2);
rhs_A=sum(sis.*rho)';
rhs_V=sum(sis.*rho.*ns)';

rhs=[rhs_f;rhs_A;rhs_V];

%===============================================================
% Reduced, vexation only system
L=[AA,AV; VA,VV];
rhs=[rhs_A;rhs_V];
sol=pinv(L)*rhs;

N=size(ff,1); B=size(fA,2); % Sizes of blocks
f=zeros(size(rho,1),1);
A=sol(1:B);
V=sol(B+1:2*B);

z_vex_sol=0;
k=prod(size(rho))-size(V,1)-size(A,1);
for n=1:size(f)
    for b=1:size(V)
    z_vex_sol=z_vex_sol+sis(n,b)*(f(n)+A(b)+V(b)*n-rho(n,b))^2/k; % Calculate z=chi^2/k  
    end
end


fprintf("====================================\n");
fprintf("Quality of fit with vexation ONLY...\n");
  z_vex_sol=z_vex_sol
  sigma=(z_vex_sol-1)*sqrt(k/2)
  rejection_confidence=( 1+erf(sigma/sqrt(2)) )/2
  one_minus_conf=erfc(sigma/sqrt(2))/2
fprintf("\n");

dA=A.*(randn(size(A))*0.001);
dV=V.*(randn(size(V))*0.001);
fprintf("Perturbation test in opposite directions...\n");
  %z_sol_pert1=resid(f*0,A+dA,V+dV,rho,sis)
  %z_sol_pert2=resid(f*0,A-dA,V-dV,rho,sis)
fprintf("\n");



% ===============================================================
%Solve for parameters
L=[ff,fA,fV; Af,AA,AV; Vf,VA,VV];
rhs=[rhs_f;rhs_A;rhs_V];
sol=pinv(L)*rhs;

N=size(ff,1); B=size(fA,2); % Sizes of blocks
f=[zeros(fnmin,1);sol(1:N)]; % Put in any zeroed out f's
A=sol(N+1:N+B);
V=sol(N+B+1:N+2*B);
figure(4); plot(f,'rx-')

z_sol=0;
k=prod(size(rho))-size(V,1)-size(f,1)-size(A,1);
for n=1:size(f)
    for b=1:size(V)
    z_sol=z_sol+sis(n,b)*(f(n)+A(b)+V(b)*n-rho(n,b))^2/k; % Calculate z=chi^2/k  
    end
end


fprintf("====================================\n");
fprintf("Quality of fit with frustration...\n");
  z_sol=z_sol
  sigma=(z_sol-1)*sqrt(k/2)
  rejection_confidence=( 1+erf(sigma/sqrt(2)) )/2
  one_minus_conf=erfc(sigma/sqrt(2))/2
fprintf("\n");

% Test some perturbation
df=f.*(randn(size(f))*0.001);
dA=A.*(randn(size(A))*0.001);
dV=V.*(randn(size(V))*0.001);
fprintf("Perturbation test in opposite directions...\n");
    
  %z_sol_pert1=resid(f+df,A+dA,V+dV,rho,sis)
  %z_sol_pert2=resid(f-df,A-dA,V-dV,rho,sis)
fprintf("\n");








