%% Polygon
pdepoly([1 0.75 -0.75 -1 -1 -0.75 0.75 1], [0.25 0.5 0.5 0.25 -0.25 -0.5 -0.5 -0.25])
% pdepoly([0 0 pi pi], [0 pi pi 0]);

%% Computation
num_triangles = size(t, 2);
num_vertices = size(p, 2);
A = zeros(num_vertices, num_vertices);
B = zeros(num_vertices, num_vertices);

for i_triangle=1:num_triangles
    indeces = t(1:3, i_triangle);
    % Add to A
    A(indeces, indeces) = A(indeces, indeces) + IntMatrix( p(:, indeces) );
    % indeces
    % A(indeces, indeces)

    % Add to B
    B(indeces, indeces) = B(indeces, indeces) + MassMatrix( p(:, indeces) );
end

[V, D] = eig(A,B); % V function, D values

% Initital conditions
c = 0.3;
g = ( 1/2 * cos(3*p(1, :)) .* cos(4*p(2,:)) )';
% g = ones(num_vertices,1);
h = ( -5*sqrt(3)*c/2 * cos(3*p(1,:)) .* cos(4*p(2,:)) )';
% h = ones(num_vertices,1);
g_hat = V \ g;
h_hat = V \ h;

dt = 0.1;
t_max = 10;
num_timesteps = 1 + (t_max / dt);
error = zeros(num_timesteps, 1);
time_vector = zeros(num_timesteps, 1);

for i=1:num_timesteps
    time = (i-1)*dt;
    lambda = diag(D);
    u_hat_t = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
    u_t = V * u_hat_t;
    u_t_exact = ( cos(3*p(1, :)) .* cos(4*p(2,:)) .* cos(5*c*time + pi/3) )';

    time_vector(i) = time;
    error(i) = max(abs(u_t - u_t_exact));
end
plot(time_vector, error);


% time = 10;
% lambda = diag(D);
% u_hat_t = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
% u_t = V * u_hat_t;
% u_t_exact = ( cos(3*p(1, :)) .* cos(4*p(2,:)) .* cos(5*c*time + pi/3) )';
% 
% pdeplot(p, [], t, 'XYData', u_t, 'ZData', u_t, 'Mesh', 'off');
% % pdeplot(p, [], t, 'XYData', u_t_exact, 'ZData', u_t_exact, 'Mesh', 'off');
% % pdeplot(p, [], t, 'XYData', g, 'ZData', g, 'Mesh', 'off');
% 
% colormap turbo

