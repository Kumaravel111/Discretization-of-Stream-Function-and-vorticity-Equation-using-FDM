clc;       % Clear the command window
clear;     % Clear all variables from the workspace
close all; % Close all figure windows

% Grid dimensions
m = 301;    % number of grid points along  x-direction
n = 301;    % number of grid points along  y-direction
dx = 1.0 / (m - 1);  % Grid size in x-direction
dy = 1.0 / (n - 1);  % Grid size in y-direction
b = dx / dy;    
re = 1000;     % Reynolds number

% matrics for storing variables
psi = zeros(m, n);    % Stream function
omega = zeros(m, n);  % Vorticity
u = zeros(m, n);      % Velocity in x-direction
v = zeros(m, n);      % Velocity in y-direction
psi_prev = zeros(m, n); 
omega_prev = zeros(m, n);
ep = 1.0;               % Error for psi
eo = 1.0;               % Error for omega
Convergence = 1e-8; % Convergence criterion
itr = 0;            % Iteration count

% Apply boundary conditions
for j = 1:n
    for i = 1:m
        if j == 1            % Bottom wall
            u(i, j) = 0.0;
            v(i, j) = 0.0;
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i, j + 1)) / dy^2;

        elseif j == n         % Top wall (moving lid)
            u(i, j) = 1.0;
            v(i, j) = 0.0;
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i, j - 1)) / dy^2 - 2 / dy;
        elseif i == 1         % Left wall
            u(i, j) = 0.0;
            v(i, j) = 0.0;
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i + 1, j)) / dx^2;
        elseif i == m         % Right wall
            u(i, j) = 0.0;
            v(i, j) = 0.0;
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i - 1, j)) / dx^2;
        end
    end
end

% using gauss seidel iterative method solve the vorticity and stream function

while (ep > Convergence || eo > Convergence)
    psi_prev = psi;             % replace previous psi and omega values with current
    omega_prev = omega;

    %  stream function 
    for j = 2:(n - 1)
        for i = 2:(m - 1)
            psi(i, j) = (0.5 / (1 + b^2)) * ...
                (psi(i + 1, j) + psi(i - 1, j) + ...
                 b^2 * (psi(i, j + 1) + psi(i, j - 1)) + ...
                 dx^2 * omega(i, j));
        end
    end

    %  vorticity

    for j = 2:(n - 1)
        for i = 2:(m - 1)
            omega(i, j) = (0.5 / (1 + b^2)) * ...
                ((1 - (psi(i, j + 1) - psi(i, j - 1)) * b * re / 4) * omega(i + 1, j) + ...
                 (1 + (psi(i, j + 1) - psi(i, j - 1)) * b * re / 4) * omega(i - 1, j) + ...
                 (1 + (psi(i + 1, j) - psi(i - 1, j)) * re / (4 * b)) * b^2 * omega(i, j + 1) + ...
                 (1 - (psi(i + 1, j) - psi(i - 1, j)) * re / (4 * b)) * b^2 * omega(i, j - 1));
        end
    end

% update boundary conditions

for j = 1:n
    for i = 1:m
        if j == 1 % Bottom wall
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i, j + 1)) / dy^2;
        elseif j == n % Top wall (moving lid)
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i, j - 1)) / dy^2 - 2 / dy;
        elseif i == 1 % Left wall
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i + 1, j)) / dx^2;
        elseif i == m % Right wall
            psi(i, j) = 0.0;
            omega(i, j) = 2 * (psi(i, j) - psi(i - 1, j)) / dx^2;
        end
    end
end

    % Compute errors

    ep = sqrt(sum((psi(:) - psi_prev(:)).^2) / (m * n));
    eo = sqrt(sum((omega(:) - omega_prev(:)).^2) / (m * n));
    itr = itr + 1;

    % Print iteration information
    fprintf('Iteration: %d, ep: %.6f, eo: %.6f\n', itr, ep, eo);
end

% Compute velocity components

for j = 2:(n - 1)
    for i = 2:(m - 1)
        u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * dy);
        v(i, j) = -(psi(i + 1, j) - psi(i - 1, j)) / (2 * dx);
    end
end

% Plotting streamlines
figure(1);
contourf(linspace(0, 1, m), linspace(0, 1, n), psi', 10);
colorbar;
caxis([min(psi(:)), max(psi(:))]);
xlabel('x');
ylabel('y');
title('Streamlines, Re = 1000');

% Plotting velocity magnitude contours
figure(2);
contourf(linspace(0, 1, m), linspace(0, 1, n), sqrt(u.^2)', 5);
colorbar;
xlabel('x');
ylabel('y');
title('Velocity u , Re = 1000');

% Plotting v-velocity along the vertical centerline
figure(3);
plot(linspace(0, 1, m), v(:, round(n / 2)), 'LineWidth', 2); % Swap axes
xlabel('x');
ylabel('v-velocity');
title('v-velocity , Re = 1000');
grid on;

% Plotting u-velocity along the vertical centerline
figure(4);
plot(u(round(m / 2), :), linspace(0, 1, n), 'LineWidth', 2); % Swap axes
xlabel('u-velocity');
ylabel('y');
title('u-velocity , Re = 1000');
grid on;

% Plot velocity vectors
figure(5);
[X, Y] = meshgrid(linspace(0, 1, m), linspace(0, 1, n));
quiver(X, Y, u', v');  % Quiver plot for velocity vectors
xlabel('x');
ylabel('y');
title('Velocity Vectors , Re = 1000');
axis equal;  
grid on;