function fnm = get_F_nm(n, m, zx, zy)
    % zx = qx * lb, zy = qy * lb, 
    z_square = zx^2 + zy^2;
    fnm = exp(-z_square / 4);
    if n >= m
        fnm = fnm * sqrt(factorial(m) / factorial(n));
        fnm = fnm * ((1j * zx - zy) / sqrt(2))^(n - m);
        fnm = fnm * laguerreL(m, n - m, z_square / 2);
        
    else
        fnm = fnm * sqrt(factorial(m) / factorial(n));
        fnm = fnm * ((1j * zx + zy) / sqrt(2))^(m - n);
        fnm = fnm * laguerreL(n, m - n, z_square / 2);
        
    end
end