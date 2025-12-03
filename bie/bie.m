%% Exercise 1

N = 300; % number of points on the boundary
tvec = linspace(-pi + 2*pi/N, pi, N); % angle vector in polar coord
rvec = 3 + cos(4*tvec+pi); % radius vector in polar coord
rprimvec = - 4*sin(4*tvec+pi); 
rbisvec = - 16*cos(4*tvec+pi);
y1 = rvec .* cos(tvec);
y2 = rvec .* sin(tvec);
nu1 = rvec .* cos(tvec) + rprimvec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimvec .* cos(tvec);
nu1 = nu1 ./ sqrt( rvec.^2 + rprimvec.^2 );
nu2 = nu2 ./ sqrt( rvec.^2 + rprimvec.^2 );

Kmat = zeros(N, N); % Kernel matrix
for i=1:N
    for j=1:N
        if i == j % diagonal case -> limit
            Kmat(i,j) = 1 / (2*pi) * ( ...
                rprimvec(i)^2 - rbisvec(i) * rvec(i) / 2 + rvec(i)^2 / 2 ...
                ) / (rprimvec(i)^2 + rvec(i)^2)^(3/2);
        else % the regular kernel
            Kmat(i,j) = 1 / (2*pi) * ( ...
                (y1(j) - y1(i)) * nu1(j) + (y2(j) - y2(i)) * nu2(j) ...
                ) / ( (y1(j) - y1(i))^2 + (y2(j) - y2(i))^2 );
        end
    end
end

% Plot kernel
% imagesc(tvec, tvec, Kmat.')
% axis xy
% colorbar

dsdt = sqrt(rprimvec.^2 + rvec.^2)'; % diagonal of Sigma matrix
% Boundary confition
gvec = ( exp( (y1 + 0.3*y2) / 3 ) .* sin( (0.3*y1 - y2) / 3 ) )';
% Solve (I/2 + 2pi/N K Sigma)h = g
hvec = ( eye(N)/2 + 2*pi/N * Kmat * diag(dsdt) ) \ gvec;

%% Compute double layer potential in D+ = solution

delta = 0.01; % delta-x = delta-y discretization
M = 8 / delta; % number of points per axis
x1field = linspace(-4, 4, M);
x2field = linspace(-4, 4, M);
ufield = zeros(M, M); % numeric solution
u0 = zeros(M, M); % analytic solution
error = -15 * ones(M, M);

for ix1 = 1:M
    for ix2 = 1:M
        x1 = x1field(ix1);
        x2 = x2field(ix2);

        % Parametrize boundary:
        t = angle( complex(x1, x2) );
        radius = 3 + cos(4*t+pi);

        if (x1^2 + x2^2) < radius^2 % Only compute within D+
            % Compute integral( del_nu Phi(y-x) h/y) ds(y) ) with trapezoid
            phivec = 1 / (2*pi) * ( (y1 - x1) .* nu1 + (y2 - x2) .* nu2 ) ./ ( (y1 - x1).^2 + (y2 - x2).^2 );
            ufield(ix1, ix2) = (phivec * (hvec .* dsdt)) * 2*pi/N;
            
            % Compute error against exact solution
            u0(ix1, ix2) = exp( (x1 + 0.3*x2) / 3 ) .* sin( (0.3*x1 - x2) / 3 );
            error(ix1, ix2) = log10( abs( u0(ix1, ix2) - ufield(ix1, ix2) ) );

        end
    end
end

% Plot

imagesc(x1field, x2field, error.');
clim([-15 0]);

% imagesc(x1field, x2field, ufield.');
% clim([-2 2]);

axis xy
xlabel("x");
ylabel("y");
colormap turbo
pbaspect([1 1 1]);
cb = colorbar;
cb.Label.String = "log10|u-u_0|";

%% Exercise 2

% Decompose indicies into even and odd
even = ( 1:1:(N/2) ) * 2;
odd = -1 + (1:1:(N/2)) * 2;
vvec = zeros(N, 1); % Define v to be calculated on boundary

% Compute v(x) at even x using odd y
for i = even
    imvec_even = 1 / (2*pi) * ( (y1(odd) - y1(i)) .* nu2(odd) - (y2(odd) - y2(i)) .* nu1(odd) ) ./ ( (y1(odd) - y1(i)).^2 + (y2(odd) - y2(i)).^2 );
    vvec(i) = (imvec_even * (hvec(odd) .* dsdt(odd))) * 2*pi/(N/2);
end

% Compute v(x) at odd x using even y
for i = odd
    imvec_odd = 1 / (2*pi) * ( (y1(even) - y1(i)) .* nu2(even) - (y2(even) - y2(i)) .* nu1(even) ) ./ ( (y1(even) - y1(i)).^2 + (y2(even) - y2(i)).^2 );
    vvec(i) = (imvec_odd * (hvec(even) .* dsdt(even))) * 2*pi/(N/2);
end

% Plot numeric v against analyic v
v_exact = ( exp( (y1 + 0.3*y2) / 3 ) .* cos( (0.3*y1 - y2) / 3 ) )';
plot(tvec, vvec, tvec, v_exact);
legend("numeric", "analytic");
xlabel("t");
ylabel("v");


%% Exercise 3

% Compute double layer potential using v(x)

% Same as in exercise 1
delta = 0.01;
M = 8 / delta;
z1field = linspace(-4, 4, M);
z2field = linspace(-4, 4, M);
gcq_ufield = zeros(M, M);
u0 = zeros(M, M);
error = -15 * ones(M, M);

for iz1 = 1:M
    for iz2 = 1:M
        z1 = z1field(iz1);
        z2 = z2field(iz2);

        % Parametrize boundary
        t = angle( complex(z1, z2) );
        radius = 3 + cos(4*t+pi);

        if (z1^2 + z2^2) < radius^2 % Only compute within D+
            
            % Helper variables
            x1 = y1 - z1;
            x2 = y2 - z2;
            norm_x = x1.^2 + x2.^2;

            % Compute Re{(i*Ip) / (i*Iq)} = Re{(Ip*conj(Iq))} / |Iq|^2
            % Where Ip = integral(p(y)dy), Iq = integral(q(y)dy)

            % p = p1 + i*p2
            p1 = ( ...
                nu1 .* (gvec' .* x1 + vvec' .* x2) + ...
                nu2 .* (gvec' .* x2 - vvec' .* x1) ...
                ) ./ norm_x;
            p2 = ( ...
                nu1 .* (vvec' .* x1 - gvec' .* x2) + ...
                nu2 .* (gvec' .* x1 + vvec' .* x2) ...
                ) ./ norm_x;
            
            % q = q1 + i*q2
            q1 = ( nu1.*x1 + nu2.*x2 ) ./ norm_x;
            q2 = ( nu2.*x1 - nu1.*x2) ./ norm_x;
            
            % Ip = Ip1 + i*Ip2
            Ip1 = (p1 * dsdt) * 2*pi/N;
            Ip2 = (p2 * dsdt) * 2*pi/N;
            
            % Iq = Iq1 + i*Iq2
            Iq1 = (q1 * dsdt) * 2*pi/N;
            Iq2 = (q2 * dsdt) * 2*pi/N;
            
            % Final computation
            gcq_ufield(iz1, iz2) = (Ip1 * Iq1 + Ip2 * Iq2) / (Iq1^2 + Iq2^2);
            
            % Compute error numeric vs analytic
            u0(iz1, iz2) = exp( (z1 + 0.3*z2) / 3 ) .* sin( (0.3*z1 - z2) / 3 );
            error(iz1, iz2) = log10( abs( u0(iz1, iz2) - gcq_ufield(iz1, iz2) ) );
        end
    end
end

% Plot
imagesc(z1field, z2field, error.');
clim([-15 0]);

% imagesc(z1field, z2field, gcq_ufield.');
% clim([-2 2]);

axis xy
xlabel("x");
ylabel("y");
colormap turbo
pbaspect([1 1 1]);
cb = colorbar;
cb.Label.String = "log10|u-u_0|";

% Check maximum errors
sorted_error = sort(error(:), 'descend');
sorted_error(1:10)
