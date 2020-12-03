function beta = instTemp(t,z,charges,domain,eVector)
	% This function uses the expression for the normal temperature for the
	% circular domain and evaluates that expression for single time points.
	% There is no time averaging used here, that is done in another function.

    n = length(charges); % Number of vortices
    tn = length(t); % Number of time points
    beta = zeros(1,tn); % Beta = 1/(k_B T)

	if domain(1) == "circle"
		% There are two conserved quantities in the circle, so we have the formula being more complicated
		% The conserved quantaties are E and L.
		for q = 1:tn
			zq = z(q,:);
			[top1, bot1, top2, bot2] = circleBeta(zq,charges);
			beta(q) = (top1/bot1) - (top2/bot2);
			beta(q) = beta(q);
		end
	elseif domain(1) == "noval"
		% There is only one conserved quantaty in the Neuman oval
		% The conserved quantaty is E
		A = str2double(domain(2)); % a paramater for the Neumann oval
		Q = str2double(domain(3)); % q paramater for the Neumann oval
		z = noval2circ(z,A,Q); % coordiante transform to the circle
		for q = 1:tn
			zq = z(q,:);
			if nargin == 5
				[top1, bot1, top2, bot2] = novalBeta(zq,charges,A,Q,eVector);
			else
				[top1, bot1, top2, bot2] = novalBeta(zq,charges,A,Q,ones(1,n));
			end
			betaq = (top1/bot1) - (top2/bot2);
			beta(q) = betaq;
		end
	end
end
