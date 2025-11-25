function error = Error(c, V, D, p, g_hat, h_hat, time)
    lambda = diag(D);
    u_hat_t = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
    u_t = V * u_hat_t;

    u_t_exact = ( cos(3*p(1, :)) .* cos(4*p(2,:)) .* cos(5*c*time + pi/3) )';
    error = max(abs(u_t - u_t_exact));

end