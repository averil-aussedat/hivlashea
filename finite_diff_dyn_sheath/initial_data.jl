const mu = 1.0/100;

""" Initial datum and boundary """
fi_0(x, v) = ones(size(x)) * exp.(-0.5 .* (v').^2)./sqrt(2pi);
fe_0(x, v) = ones(size(x)) * sqrt(mu).*exp.(-0.5 * mu .* (v').^2)./sqrt(2*pi);
E0(x) = 0.0*x;

