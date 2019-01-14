function bs = dsphbesu(n, x, type)
    % retuns the derivate of spherical Bessel function
    bs = (n.*sphbesu(n-1,x,type) - (n+1).*sphbesu(n+1,x,type))./(2.*n+1);
end

