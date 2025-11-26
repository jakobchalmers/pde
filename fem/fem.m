%% Polygon
pdepoly([1 0.75 -0.75 -1 -1 -0.75 0.75 1], [0.25 0.5 0.5 0.25 -0.25 -0.5 -0.5 -0.25])
% pdepoly([0 0 pi pi], [0 pi pi 0]);

%% Exercise 1
% save("p_h4", "p");

num_triangles = size(t, 2);
num_vertices = size(p, 2);
A = zeros(num_vertices, num_vertices);
B = zeros(num_vertices, num_vertices);

for i_triangle=1:num_triangles
    indeces = t(1:3, i_triangle); % Indices for nodes in the triangle

    % Add contributions from the triangle to A
    A(indeces, indeces) = A(indeces, indeces) + StiffnessMatrix( p(:, indeces) );

    % Add contributions from the triangle to B
    B(indeces, indeces) = B(indeces, indeces) + MassMatrix( p(:, indeces) );
end

[V, D] = eig(A,B); % Solve Rayleigh Ritz
lambda = diag(D); % vector of eigenvalues
% save('V_h4.mat', 'V');
% save('D_h4.mat', 'D');

%%

% Initital conditions
c = 0.3;
g = ( 1/2 * cos(3*p(1, :)) .* cos(4*p(2,:)) )';
h = ( -5*sqrt(3)*c/2 * cos(3*p(1,:)) .* cos(4*p(2,:)) )';
g_hat = V \ g;
h_hat = V \ h;
% save("g_hat_h4.mat", "g_hat");
% save("h_hat_h4.mat", "h_hat");

% Solve and plot u
time = 10;
u_hat_t = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
u_t = V * u_hat_t;
% u_t_exact = ( cos(3*p(1, :)) .* cos(4*p(2,:)) .* cos(5*c*time + pi/3) )';
pdeplot(p, [], t, 'XYData', u_t, 'ZData', u_t, 'Mesh', 'off');
colormap turbo


%% Exercise 2
c = 0.3;
% Discretize time
dt = 1;
t_max = 10;
num_timesteps = 1 + (t_max / dt);
time_vector = zeros(num_timesteps, 1);

% Error vector for each triangulation
error_h1 = zeros(num_timesteps, 1);
error_h2 = zeros(num_timesteps, 1);
error_h3 = zeros(num_timesteps, 1);
error_h4 = zeros(num_timesteps, 1);

for i=1:num_timesteps
    time = (i-1)*dt;
    time_vector(i) = time;
    
    % Compute maximum absolute error for each mesh
    error_h1(i) = Error(c, load("V_h1.mat").V, load("D_h1.mat").D, load("p_h1.mat").p, load("g_hat_h1.mat").g_hat, load("h_hat_h1.mat").h_hat, time);
    error_h2(i) = Error(c, load("V_h2.mat").V, load("D_h2.mat").D, load("p_h2.mat").p, load("g_hat_h2.mat").g_hat, load("h_hat_h2.mat").h_hat, time);
    error_h3(i) = Error(c, load("V_h3.mat").V, load("D_h3.mat").D, load("p_h3.mat").p, load("g_hat_h3.mat").g_hat, load("h_hat_h3.mat").h_hat, time);
    error_h4(i) = Error(c, load("V_h4.mat").V, load("D_h4.mat").D, load("p_h4.mat").p, load("g_hat_h4.mat").g_hat, load("h_hat_h4.mat").h_hat, time);

end
% Plot Error vs. time for each mesh
plot(time_vector, error_h1, time_vector, error_h2, time_vector, error_h3, time_vector, error_h4);
xlabel("t")
ylabel("Error")
title("Error vs. time")
legend('h1', 'h2', 'h3', 'h4')

% Compute mean errors over time for each mesh
mean_error = [
    mean(error_h1);
    mean(error_h2);
    mean(error_h3);
    mean(error_h4);
    ];
save("mean_error.mat", "mean_error");

%% 
% Compute approximate triangle side length h
num_nodes = [175; 655; 2533; 9961];
h = sqrt(2 * pi^2 ./ num_nodes);

% Find N in Ordo(h^N)
mean_error = load("data/mean_error.mat").mean_error;
fit = polyfit(log(h), log(mean_error), 1); % Fit a line
N = fit(2) % slope of the line


%% Exercise 3 a)
% Compute the index of the node closest to (1, 0)
p0 = [1; 0]; 
distance = (p(1,:) - p0(1)).^2 + (p(2,:) - p0(2)).^2;
[dist, node_index] = min(distance);
dist
node_index

%%

% Initital conditions
c = 0.3;
g = exp( -100*( (p(1,:) - 0.5).^2 + p(2, :).^2 ) )';
h = p(1,:)' * 0;
g_hat = V \ g;
h_hat = V \ h;

%%

% Solve
dt = 0.01;
t_max = 2;
num_timesteps = 1 + (t_max / dt);
amp_vector = zeros(num_timesteps, 1); % u_t(x=1, y=0) over time

% Find u_t(x=1, y=0) for each time
for i=1:num_timesteps
    time = (i-1)*dt;
    u_hat_t = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
    u_t = V * u_hat_t;
    
    amp_vector(i) = u_t(node_index);
end

% Find time T of hitting right edge and deduced propagation speed
[max_amp, time_index] = max(amp_vector);
T = (time_index-1)*dt
c_estimate = 0.5 / c

pdeplot(p, [], t, 'XYData', u_t, 'ZData', u_t, 'Mesh', 'off');
colormap turbo

%% 3 b)

% Compute u_10
time = 10;
u_hat_10 = cos(c*sqrt(lambda)*time) .* g_hat + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* h_hat;
u_10 = V * u_hat_10; % Save as IC 
del_t_u_hat_10 = - c*sqrt(lambda) .* sin(c*sqrt(lambda)*time) .* g_hat + cos(c*sqrt(lambda)*time) .* h_hat;
del_t_u_10 = V * del_t_u_hat_10; % Save as IC

% Compute v_{-9} using u_10 and del_t_u_10
time = -9;
v_hat_neg9 = cos(c*sqrt(lambda)*time) .* u_hat_10 + ( sin(c*sqrt(lambda)*time) ./ (c*sqrt(lambda)) ) .* del_t_u_hat_10;
v_neg9 = V * v_hat_neg9;


%% u_10 plot
pdeplot(p, [], t, 'XYData', u_10, 'ZData', u_10, 'Mesh', 'off');
colormap turbo

%% del_t(u_10) plot
pdeplot(p, [], t, 'XYData', del_t_u_10, 'ZData', del_t_u_10, 'Mesh', 'off');
colormap turbo

%% v_{-9} plot
pdeplot(p, [], t, 'XYData', v_neg9, 'ZData', v_neg9, 'Mesh', 'off');
colormap turbo




    
 