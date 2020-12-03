function dzdt = PointVortexPlane(t,zeta,kappa,gamma0,N,rho0,R) 

Rt = @(t) R*sqrt(1 + 2*pi*rho0*gamma0*t);
dRdt = @(t) pi*R*rho0*gamma0/sqrt(1+2*pi*rho0*gamma0*t);
Omegat = @(t) N/Rt(t)^2;

tt = @(t) (exp(R^2*rho0*gamma0*2*pi*t) - 1)/rho0/gamma0/2/pi;

dzdt = zeros(N,1);
for ii = 1:N
    realvortices = 1./(zeta(ii) - zeta);
    realvortices(ii) = 0;
    %imagevortices = conj(zeta)./(1- zeta(ii)*conj(zeta));
    dzdt(ii) = -1i*kappa*(realvortices);
end


% normal dynamics in scaled coordinates -- space scaled only
%  dzdt = (conj(dzdt) - 1i*gamma0*kappa.'.*conj(dzdt));
%  dzdt =  (dzdt - (1i*N + dRdt(t)*Rt(t))*zeta) /Rt(t)^2;

  %"Accelerated" dynamics --- unsure if this is correct yet
  dzdt = (conj(dzdt) - 1i*gamma0*kappa.'.*conj(dzdt));
  dzdt =  (dzdt - (1i*N + dRdt((t))*Rt((t)))*zeta);


% normal dynamics -- no scaling
 %dzdt = (conj(dzdt) - 1i*gamma0*kappa.'.*conj(dzdt));