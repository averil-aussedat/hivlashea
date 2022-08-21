const mu = 0.5 #1.0/100;

""" Initial datum and boundary """
#fi_0(x, v) = ones(size(x)) * exp.(-0.5 .* (v').^2)./sqrt(2*pi);
#fe_0(x, v) = ones(size(x)) * sqrt(mu).*exp.(-0.5 * mu .* (v').^2)./sqrt(2*pi);
mask(x) = 0.5*(tanh.((x+0.1)./0.1)-tanh.((x-0.1)./0.1))
#fi_0(x, v) = mask.(x).* exp.(-0.5 .* (v').^2)./sqrt(2*pi)
#fe_0(x, v) = mask.(x).* sqrt(mu).*exp.(-0.5 * mu .* (v').^2)./sqrt(2*pi)
fi_0(x, v) = mask.(x).* exp.(-0.5 .* (v').^2)./sqrt(2*pi)
fe_0(x, v) = sqrt(mu)*mask.(x).* exp.(-0.5 * mu .* (v').^2)./sqrt(2*pi);
E0(x) = 0.0*x

