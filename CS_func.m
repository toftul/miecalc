function CS=CS_func(r_in)
%----------Input  Parameters------------------------------------------------
n = 3; % Multipole order 
r =r_in;% 120/2*1e-9; %[m] Sphere radius 
S0=1%4*pi*r^2; % [m^2] geometrical cross section
I = 1.3e-3; %[W] Incident wave intensity 
lambda = 630e-9; %[m] Wavelength
f = 3e8./lambda; %[Hz] frequency
w=1.24e-6./lambda; % [eV] energy
mu = 1; % magnetic permiability of medium
mu1 = 1; % magnetic permiability of sphere
N = 1; % refractive index of medium
k = 2*pi*N./lambda; % [1/m] wavevector 
%-------------------Material Parameters------------------------------------
Material='Si';
epsi=EpsMat(w,Material);
N1 = sqrt(epsi);
%----------Calculating cross sections--------------------------------------
%Scattering cross section
Csca = 0;
CscaE = 0;
CscaM = 0;
for i=1:1:n
    Csca = Csca + 2*pi./k.^2 .* (2.*i+1) .* (abs(An(i,mu,mu1,N,N1,r,lambda)).^2 ...
    + abs(Bn(i,mu,mu1,N,N1,r,lambda)).^2);
    CscaE = CscaE + 2*pi./k.^2 .* (2.*i+1) .* (abs(An(i,mu,mu1,N,N1,r,lambda)).^2);
    CscaM = CscaM + 2*pi./k.^2 .* (2.*i+1) .* (abs(Bn(i,mu,mu1,N,N1,r,lambda)).^2);
end

%Extinction cross section
Cext = 0;
for i=1:1:n
    Cext = Cext + 2*pi./k.^2 .* (2.*i+1) .* real(An(i,mu,mu1,N,N1,r,lambda)...
        + Bn(i,mu,mu1,N,N1,r,lambda));
end
CS=Cext;
return

%-------------Plotting results--------------------------------------------

% figure(1);
% hold on
% plot(lambda*1e9,Csca/pi/r^2,'k','linewidth',2);
% plot(lambda*1e9,CscaE/pi/r^2,'r','linewidth',2);
% plot(lambda*1e9,CscaM/pi/r^2,'b','linewidth',2);
% %plot(lambda*1e9,Cext/S0,'r','linewidth',2);
% xlabel('Wavelength, nm,','Fontsize',20);
% ylabel('Scattering cross section, \mum^2 ','Fontsize',20);
% legend('Scattering ','Extinction ')
% set(gca,'fontsize',20)
