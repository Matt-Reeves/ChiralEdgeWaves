function dzdt = getVelocity(zeta,kappa,N) 


dzdt = zeros(N,1);
for ii = 1:N
    realvortices = 1./(zeta(ii) - zeta);
    realvortices(ii) = 0;
    dzdt(ii) = -1i*kappa*(realvortices);
end

dzdt = conj(dzdt);