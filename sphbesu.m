function bs = sphbesu(nu, x, type)
% returns the spherical Bessel functions jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix

[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);

switch type
    case 'j'
         bs = sqrt(pi ./(2* xm)) .* besselj(nu + 0.5, x);
    case 'y'
        bs = sqrt(pi ./(2* xm)) .* bessely(nu + 0.5, x);
    case 'h'
        bs = sqrt(pi ./(2* xm)) .* besselh(nu + 0.5, x);
end


