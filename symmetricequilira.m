function dydt = symmetricequilira(t,y,f,g)

if t == 0
dydt = [0; f*exp(y(1))-g];
else
dydt = [y(2); -y(2)/t+f*exp(y(1))-g];
end