function [top1, bot1, top2, bot2] = circleBeta(z,charges)

	n = length(charges);

	LF = -pi;
	HF = pi;

	% L partial derivatives
	Lz = LF*charges.*conj(z);
	Lzbar = LF*charges.*z;
	Lzz = zeros(n,n);
	Lzbarzbar = zeros(n,n);
	Lzbarz = zeros(n,n);
	for jj = 1:n
		Lzbarz(jj,jj) = LF*charges(jj);
	end

	% H partial derivatives
	Hz = zeros(1,n);
	Hzbar = zeros(1,n);
	Hzz = zeros(n,n);
	Hzbarzbar = zeros(n,n);
	Hzbarz = zeros(n,n);
	for ii = 1:n
		zi = z(ii);
		Hzi = -charges(ii).*charges./(zi-z);
		Hzi(ii) = 0;
		Hzi = Hzi - charges(ii).*charges.*conj(z)./(1-zi*conj(z));
		Hz(ii) = HF*sum(Hzi);

		Hzbar(ii) = conj(Hz(ii));

		for jj = 1:n
			zj = z(jj);
			if jj == ii
				Hzzii = charges(ii).*charges./(zi-z).^2;
				Hzzii(ii) = 0;
				Hzzii = Hzzii - charges(ii).*charges.*(conj(z).^2)./(1-conj(z).*zi).^2;
				Hzz(ii,ii) = HF*sum(Hzzii);

				Hzbarzbar(ii,ii) = conj(Hzz(ii,ii));

				Hzbarz(ii,ii) = -HF*charges(ii)^2 ./(1-zi.*conj(zi))^2;
			else
				Hzz(ii,jj) = -HF*charges(ii)*charges(jj)/(zi-zj)^2;

				Hzbarzbar(ii,jj) = conj(Hzz(ii,jj));

				Hzbarz(ii,jj) = -HF*charges(ii)*charges(jj)/(1-zj*conj(zi))^2;
			end
		end
	end

	del2H = 4*ones(1,n)*diag(Hzbarz);
	del2L = 4*ones(1,n)*diag(Lzbarz);
	delHdelH = 4*Hz*transpose(Hzbar);
	delLdelL = 4*Lz*transpose(Lzbar);
	delHdelL = 2*(Hz*transpose(Lzbar) + Lz*transpose(Hzbar));

	HHz = 4*(Hzbar*Hzz + Hz*Hzbarz);
	HHzbar = 4*(Hzbar*transpose(Hzbarz) + Hz*Hzbarzbar);
	LLz = 4*(Lzbar*Lzz + Lz*Lzbarz);
	LLzbar = 4*(Lzbar*transpose(Lzbarz) + Lz*Lzbarzbar);
	HLz = 2*(Lzbar*Hzz + Hz*Lzbarz + Lz*Hzbarz + Hzbar*Lzz);
	HLzbar = 2*(Lz*Hzbarzbar + Hz*Lzbarzbar + Lzbar*transpose(Hzbarz) + Hzbar*transpose(Lzbarz));

	lambda = delHdelL/delLdelL;

	lambz = (delLdelL*HLz - delHdelL*LLz)./delLdelL^2;
	lambzbar = (delLdelL*HLzbar - delHdelL*LLzbar)./delLdelL^2;

	dellambdelL = 2*(Lz*transpose(lambzbar)+Lzbar*transpose(lambz));

	% We are going to express this dot product as:
	% A o B + C o D;
	A = HHz - delHdelL*lambz - lambda*HLz;
	B = transpose(Hzbar) - lambda*(transpose(Lzbar));
	C = HHzbar - delHdelL*lambzbar - lambda*HLzbar;
	D = transpose(Hz) - lambda*(transpose(Lz));

	bigdotprod = 2*(A*B+C*D);

	top1 = del2H - lambda*del2L - dellambdelL;
	bot1 = delHdelH - lambda*delHdelL;
	top2 = bigdotprod;
	bot2 = bot1^2;

end
