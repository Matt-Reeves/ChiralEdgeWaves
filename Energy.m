function H = Energy(zeta,N)

%Calculates the energy under conformal mapping to the unit circle.
   
Hint = 0;
for jj = 1:N
    fij = abs( (zeta(jj) - zeta));
    fij = log(fij);
    fij(jj) = 0;
    Hint = Hint - sum(fij);
end
H = Hint;

end


