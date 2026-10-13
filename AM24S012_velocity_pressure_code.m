% MATLAB Code for 2D Channel Flow with Excel Input Parameters
clear all; close all; clc;

%% Read Parameters from Excel with preserved headers
params = readtable('AM24S012_inputfile.xlsx', ...
    'VariableNamingRule','preserve', ...
    'Range','A1:C10');

% Extract parameters using cleaned header name
param_col = strtrim(params.("Problem Parameters"));
L = params.Value(strcmp(param_col, "Channel length"));
H = params.Value(strcmp(param_col, "Channel height")); 
u_in = params.Value(strcmp(param_col, "Inlet velocity"));
Re = params.Value(strcmp(param_col, "Reynolds number"));
rho = params.Value(strcmp(param_col, "Density"));
Nx = params.Value(strcmp(param_col, "Number of cells in x-direction"));
Ny = params.Value(strcmp(param_col, "Number of cells in y-direction"));
dt = params.Value(strcmp(param_col, "Time step"));
max_iter = params.Value(strcmp(param_col, "Maximum number of iterations"));
tol = params.Value(strcmp(param_col, "Convergence tolerance"));
n_correctors = 2;    % Number of PISO corrector steps
p_iter = 50;         % Pressure correction iterations
omega = 1.5;         % SOR relaxation factor

%% Derived Parameters
mu = rho*u_in*H/Re;  % Dynamic viscosity
nu = mu/rho;         % Kinematic viscosity

%% Rest of the Original Code Remains Unchanged
%% Grid Generation
dx = L/Nx;
dy = H/Ny;

% Cell centers
xc = linspace(dx/2, L-dx/2, Nx);
yc = linspace(dy/2, H-dy/2, Ny);
[X, Y] = meshgrid(xc, yc);
X = X';
Y = Y';

%% Initialize Fields
u = zeros(Nx, Ny);
v = zeros(Nx, Ny);
p = zeros(Nx, Ny);
u_f = zeros(Nx+1, Ny);
v_f = zeros(Nx, Ny+1);

% Set inlet boundary condition
u_f(1,:) = u_in;

%% Main Solution Loop
for iter = 1:max_iter
    % Store old values for convergence check
    u_old = u;
    v_old = v;
    p_old = p;
    
    %% PISO Algorithm Steps
    %% Step 1: Momentum Predictor
    % Calculate H operator for u-momentum
    H_u = zeros(Nx, Ny);
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Calculate convective terms
            dudx = (u(i+1,j) - u(i-1,j))/(2*dx);
            dudy = (u(i,j+1) - u(i,j-1))/(2*dy);
            
            % Calculate diffusive terms
            d2udx2 = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/(dx^2);
            d2udy2 = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/(dy^2);
            
            % H operator (all terms except pressure gradient)
            H_u(i,j) = u(i,j) + dt*(-u(i,j)*dudx - v(i,j)*dudy + nu*(d2udx2 + d2udy2));
        end
    end
    
    % Calculate H operator for v-momentum
    H_v = zeros(Nx, Ny);
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Calculate convective terms
            dvdx = (v(i+1,j) - v(i-1,j))/(2*dx);
            dvdy = (v(i,j+1) - v(i,j-1))/(2*dy);
            
            % Calculate diffusive terms
            d2vdx2 = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/(dx^2);
            d2vdy2 = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/(dy^2);
            
            % H operator (all terms except pressure gradient)
            H_v(i,j) = v(i,j) + dt*(-u(i,j)*dvdx - v(i,j)*dvdy + nu*(d2vdx2 + d2vdy2));
        end
    end
    
    % PISO correction loop
    for corr = 1:n_correctors
        %% Step 2: Calculate face velocities
        % x-face velocities
        for i = 2:Nx
            for j = 2:Ny-1
                % Interpolate H to face and add pressure gradient
                H_face = 0.5*(H_u(i-1,j) + H_u(i,j));
                dp_dx = (p(i,j) - p(i-1,j))/dx;
                u_f(i,j) = H_face - dt*dp_dx/rho;
            end
        end
        
        % y-face velocities
        for i = 2:Nx-1
            for j = 2:Ny
                % Interpolate H to face and add pressure gradient
                H_face = 0.5*(H_v(i,j-1) + H_v(i,j));
                dp_dy = (p(i,j) - p(i,j-1))/dy;
                v_f(i,j) = H_face - dt*dp_dy/rho;
            end
        end
        
        %% Step 3: Calculate divergence and solve pressure correction
        % Calculate divergence of velocity field
        div = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
                div(i,j) = (u_f(i+1,j) - u_f(i,j))/dx + (v_f(i,j+1) - v_f(i,j))/dy;
            end
        end
        
        % Solve pressure correction equation
        p_corr = zeros(Nx, Ny);
        for k = 1:p_iter
            p_old_iter = p_corr;
            
            for i = 2:Nx-1
                for j = 2:Ny-1
                    % Coefficients for pressure equation
                    aE = dt/(rho*dx^2);
                    aW = dt/(rho*dx^2);
                    aN = dt/(rho*dy^2);
                    aS = dt/(rho*dy^2);
                    aP = 2*dt/(rho*dx^2) + 2*dt/(rho*dy^2);
                    
                    % Source term (negative divergence)
                    b = -div(i,j);
                    
                    % SOR iteration
                    p_corr(i,j) = (1-omega)*p_corr(i,j) + omega/aP * ...
                                 (aE*p_corr(i+1,j) + aW*p_corr(i-1,j) + ...
                                  aN*p_corr(i,j+1) + aS*p_corr(i,j-1) + b);
                end
            end
            
            % Pressure correction boundary conditions
            p_corr(1,:) = p_corr(2,:);        % Zero gradient at inlet
            p_corr(Nx,:) = 0;                 % Zero pressure at outlet
            p_corr(:,1) = p_corr(:,2);        % Zero gradient at bottom wall
            p_corr(:,Ny) = p_corr(:,Ny-1);    % Zero gradient at top wall
            
            % Check convergence of pressure iteration
            if max(max(abs(p_corr - p_old_iter))) < 1e-4
                break;
            end
        end
        
        %% Step 4: Correct pressure
        p = p + p_corr;
        
        %% Step 5: Correct face velocities
        % Correct x-face velocities
        for i = 2:Nx
            for j = 2:Ny-1
                dp_corr_dx = (p_corr(i,j) - p_corr(i-1,j))/dx;
                u_f(i,j) = u_f(i,j) - dt*dp_corr_dx/rho;
            end
        end
        
        % Correct y-face velocities
        for i = 2:Nx-1
            for j = 2:Ny
                dp_corr_dy = (p_corr(i,j) - p_corr(i,j-1))/dy;
                v_f(i,j) = v_f(i,j) - dt*dp_corr_dy/rho;
            end
        end
        
        %% Step 6: Correct cell-center velocities
        % Correct x-velocity
        for i = 2:Nx-1
            for j = 2:Ny-1
                dp_corr_dx = (p_corr(i+1,j) - p_corr(i-1,j))/(2*dx);
                u(i,j) = H_u(i,j) - dt*dp_corr_dx/rho;
            end
        end
        
        % Correct y-velocity
        for i = 2:Nx-1
            for j = 2:Ny-1
                dp_corr_dy = (p_corr(i,j+1) - p_corr(i,j-1))/(2*dy);
                v(i,j) = H_v(i,j) - dt*dp_corr_dy/rho;
            end
        end
    end
    
    %% Apply boundary conditions
    % Inlet (x = 0)
    u(1,:) = 2*u_in - u(2,:);  % Extrapolation to maintain u_f(1,:) = u_in
    v(1,:) = 0;
    
    % Outlet (x = L)
    u(Nx,:) = u(Nx-1,:);       % Zero gradient (∂u/∂x = 0)
    v(Nx,:) = v(Nx-1,:);       % Zero gradient (∂v/∂x = 0)
    
    % Bottom wall (y = 0)
    u(:,1) = 0;                % No-slip
    v(:,1) = 0;                % No-slip
    
    % Top wall (y = H)
    u(:,Ny) = 0;               % No-slip
    v(:,Ny) = 0;               % No-slip
    
    % Face velocity boundary conditions
    u_f(1,:) = u_in;           % Inlet
    u_f(Nx+1,:) = u_f(Nx,:);   % Outlet
    v_f(:,1) = 0;              % Bottom wall
    v_f(:,Ny+1) = 0;           % Top wall
    
    %% Check convergence
    u_error = max(max(abs(u - u_old)));
    v_error = max(max(abs(v - v_old)));
    p_error = max(max(abs(p - p_old)));
    
    max_error = max([u_error, v_error, p_error]);
    
    if mod(iter, 100) == 0
        fprintf('Iteration %d: max error = %.6e\n', iter, max_error);
    end
    
    if max_error < tol
        fprintf('Solution converged after %d iterations\n', iter);
        break;
    end
end

%% Post-processing
% Calculate velocity magnitude
vel_mag = sqrt(u.^2 + v.^2);

% Plot velocity magnitude contour
figure(1);
contourf(X, Y, vel_mag, 20);
colormap(jet);
colorbar;
title('Velocity Magnitude');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Plot pressure contour
figure(2);
contourf(X, Y, p, 20);
colormap(jet);
colorbar;
title('Pressure Distribution');
xlabel('x-coordinate');
ylabel('y-coordinate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% u velocity profile

% Plot u-velocity profiles at different x-positions
figure(6);
set(gcf, 'Position', [100 100 1200 600]);

% Define the x-positions to plot (in units)
x_positions = 1:10;

% Create subplots - 2 rows, 5 columns
for i = 1:length(x_positions)
    subplot(2, 5, i);
    hold on;
    
    x_pos = x_positions(i);
    
    % Find the closest index in xc to the desired x-position
    [~, idx] = min(abs(xc - x_pos));
    
    % Plot u-velocity vs. y at this x-position
    plot(u(idx,:), yc, 'b-o', 'LineWidth', 1.5);
    
    % Plot analytical solution for comparison
    y_analytical = linspace(0, H, 100);
    u_analytical = 6*u_in*y_analytical.*(H - y_analytical)/H^2;
    plot(u_analytical, y_analytical, 'r--', 'LineWidth', 1.5);
    
    grid on;
    title(sprintf('x = %d units', x_pos));
    xlabel('u-velocity');
    ylabel('y-coordinate');
    
    if i == 1
        legend('Numerical', 'Analytical', 'Location', 'best');
    end
    
    % Set consistent axis limits for better comparison
    xlim([0 max(u_analytical)*1.1]);
    ylim([0 H]);
end

sgtitle('U-velocity Profiles at Different x-positions', 'FontSize', 14);

% Alternatively, plot all profiles on a single graph
figure(7);
hold on;

% Create a colormap for different x-positions
colors = jet(length(x_positions));
legend_entries = cell(length(x_positions), 1);

% Plot each profile with a different color
for i = 1:length(x_positions)
    x_pos = x_positions(i);
    [~, idx] = min(abs(xc - x_pos));
    plot(u(idx,:), yc, 'LineWidth', 1.5, 'Color', colors(i,:));
    legend_entries{i} = sprintf('x = %d', x_pos);
end

% Add analytical solution
y_analytical = linspace(0, H, 100);
u_analytical = 6*u_in*y_analytical.*(H - y_analytical)/H^2;
plot(u_analytical, y_analytical, 'k--', 'LineWidth', 2);
legend_entries{end+1} = 'Analytical (fully developed)';

grid on;
title('Comparison of U-velocity Profiles');
xlabel('u-velocity');
ylabel('y-coordinate');
legend(legend_entries, 'Location', 'eastoutside');
% % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%   Entry length

%% Calculate Entry Length (L_E) where flow becomes fully developed
% A flow is considered fully developed when velocity profile stops changing 
% in the streamwise direction (du/dx ≈ 0)

% Calculate du/dx using the gradient function
[dudx, ~] = gradient(u, dx, dy);

% Calculate average absolute du/dx at each x-position
du_dx_avg = zeros(Nx, 1);
for i = 1:Nx
    du_dx_avg(i) = mean(abs(dudx(i,:)));
end

% Plot the rate of velocity profile change along the channel
figure(6);
plot(xc, du_dx_avg, 'b-o', 'LineWidth', 1.5);
grid on;
title('Rate of Change of u-velocity Along Channel (du/dx)');
xlabel('x-coordinate');
ylabel('Average |du/dx|');

% Define threshold for "fully developed" flow
% Typically a small percentage of maximum gradient
threshold = 0.01 * max(du_dx_avg);  % 1% of maximum rate of change

% Find position where du/dx falls below threshold
entry_indices = find(du_dx_avg < threshold);
if ~isempty(entry_indices)
    entry_index = entry_indices(1);  % First position meeting criteria
    entry_length = xc(entry_index);
    
    % Mark the entry length on the plot
    hold on;
    plot([entry_length, entry_length], [0, du_dx_avg(entry_index)], 'r--', 'LineWidth', 1.5);
    plot(entry_length, du_dx_avg(entry_index), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(entry_length+0.1, 2*du_dx_avg(entry_index), ['L_E = ', num2str(entry_length, '%.2f')], 'FontSize', 12);
    
    % Display results
    fprintf('Entry length (L_E) = %.4f units\n', entry_length);
    fprintf('Entry length ratio (L_E/H) = %.4f\n', entry_length/H);
    
    % Compare with theoretical prediction for laminar flow
    theoretical_Le = 0.06 * Re * H;  % From laminar flow theory
    fprintf('Theoretical entry length (0.06 * Re * H) = %.4f units\n', theoretical_Le);
else
    fprintf('Flow not fully developed within simulation domain\n');
end

% Create visualization showing profile development
figure(7);
colormap jet;
subplot(2,1,1);
contourf(X, Y, dudx, 20);
colorbar;
hold on;
if exist('entry_length', 'var')
    plot([entry_length, entry_length], [0, H], 'r--', 'LineWidth', 2);
    text(entry_length+0.1, H/2, 'L_E', 'Color', 'r', 'FontSize', 12);
end
title('Spatial Derivative of u-velocity (du/dx)');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Compare velocity profiles at different x-positions
subplot(2,1,2);
hold on;
x_positions = [1, round(entry_index/3), round(2*entry_index/3), entry_index, Nx];
legends = cell(length(x_positions)+1, 1);
colors = jet(length(x_positions));

for i = 1:length(x_positions)
    pos = x_positions(i);
    if pos > Nx, pos = Nx; end
    plot(u(pos,:), yc, 'LineWidth', 1.5, 'Color', colors(i,:));
    legends{i} = ['x = ', num2str(xc(pos), '%.2f')];
end

% Add analytical fully developed profile
y_analytical = linspace(0, H, 100);
u_analytical = 6*u_in*y_analytical.*(H - y_analytical)/H^2;
plot(u_analytical, y_analytical, 'k--', 'LineWidth', 2);
legends{end} = 'Analytical';

grid on;
title('U-velocity Profiles at Different x-positions');
xlabel('u-velocity');
ylabel('y-coordinate');
legend(legends, 'Location', 'best');

% % Calculate maximum velocity at each x-position to show development
% u_max = max(u, [], 2);
% figure(8);
% plot(xc, u_max, 'b-', 'LineWidth', 1.5);
% hold on;
% if exist('entry_length', 'var')
%     plot([entry_length, entry_length], [min(u_max), max(u_max)], 'r--', 'LineWidth', 1.5);
%     text(entry_length+0.1, 0.95*max(u_max), 'L_E', 'Color', 'r', 'FontSize', 12);
% end
% grid on;
% title('Maximum u-velocity Along Channel Length');
% xlabel('x-coordinate');
% ylabel('Maximum u-velocity');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  comparison with analytical ans

%% 
%% Define analytical solutions for comparison
% Calculate analytical velocity and pressure fields
u_analytical_field = zeros(size(X));
p_analytical_field = zeros(size(X));

% Calculate theoretical entry length
Le_theoretical = 0.06 * Re * H;

% Set up analytical solutions
for i = 1:Nx
    for j = 1:Ny
        % Position coordinates
        x_pos = X(i,j);
        y_pos = Y(i,j);
        
        % Fully developed velocity profile
        if x_pos >= Le_theoretical
            % Parabolic profile
            u_analytical_field(i,j) = 6 * u_in * y_pos * (H - y_pos) / H^2;
        else
            % Simple blending in entry region
            blend_factor = x_pos / Le_theoretical;
            u_fd = 6 * u_in * y_pos * (H - y_pos) / H^2;
            u_analytical_field(i,j) = u_in * (1 - blend_factor) + u_fd * blend_factor;
        end
        
        % Analytical pressure (linear decrease along channel)
        dpdx = -12 * mu * u_in / (H^2);
        p_analytical_field(i,j) = dpdx * x_pos;
    end
end

% Calculate analytical velocity magnitude
vel_mag_analytical = abs(u_analytical_field);

% Normalize pressure fields for comparison
p = p - p(1, round(Ny/2));  % Set inlet centerline pressure to zero
p_analytical_field = p_analytical_field - p_analytical_field(1, round(Ny/2));

%% Compare velocity magnitude
figure(9);
set(gcf, 'Position', [100 100 1200 500]);

% Numerical solution
subplot(1,2,1);
contourf(X, Y, vel_mag, 20);
colormap(jet);
colorbar;
title('Numerical Velocity Magnitude');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Analytical solution
subplot(1,2,2);
contourf(X, Y, vel_mag_analytical, 20);
colormap(jet);
colorbar;
title('Analytical Velocity Magnitude');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Calculate and display difference
vel_mag_diff = vel_mag - vel_mag_analytical;
max_vel_diff = max(max(abs(vel_mag_diff)));
avg_vel_diff = mean(mean(abs(vel_mag_diff)));

% figure(10);
% contourf(X, Y, vel_mag_diff, 20);
% colormap(jet);
% colorbar;
% title(['Velocity Magnitude Difference (Max: ', num2str(max_vel_diff, '%.5f'), ', Avg: ', num2str(avg_vel_diff, '%.5f'), ')']);
% xlabel('x-coordinate');
% ylabel('y-coordinate');
%% Compare pressure distribution
figure(11);
set(gcf, 'Position', [100 100 1200 500]);

% Numerical solution
subplot(1,2,1);
contourf(X, Y, p, 20);
colormap(jet);
colorbar;
title('Numerical Pressure Distribution');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Analytical solution
subplot(1,2,2);
contourf(X, Y, p_analytical_field, 20);
colormap(jet);
colorbar;
title('Analytical Pressure Distribution');
xlabel('x-coordinate');
ylabel('y-coordinate');

% Calculate and display difference
p_diff = p - p_analytical_field;
max_p_diff = max(max(abs(p_diff)));
avg_p_diff = mean(mean(abs(p_diff)));

% figure(12);
% contourf(X, Y, p_diff, 20);
% colormap(jet);
% colorbar;
% title(['Pressure Difference (Max: ', num2str(max_p_diff, '%.5f'), ', Avg: ', num2str(avg_p_diff, '%.5f'), ')']);
% xlabel('x-coordinate');
% ylabel('y-coordinate');

%% Enhanced U-velocity profile comparison
figure(13);
set(gcf, 'Position', [100 100 1200 800]);

% Define sampling positions along channel
x_positions = linspace(0.5, L-0.5, 10);
colors = jet(length(x_positions));

% 1. Plot u-velocity profiles with analytical comparison
subplot(2,2,1);
hold on;
grid on;

for i = 1:length(x_positions)
    x_pos = x_positions(i);
    [~, idx] = min(abs(xc - x_pos));
    plot(u(idx,:), yc, 'Color', colors(i,:), 'LineWidth', 1.5);
end

% Add analytical profile
y_analytical = linspace(0, H, 100);
u_analytical = 6*u_in*y_analytical.*(H - y_analytical)/H^2;
plot(u_analytical, y_analytical, 'k--', 'LineWidth', 2);

title('U-velocity Profiles at Different Locations');
xlabel('u-velocity');
ylabel('y-coordinate');
legend([cellstr(num2str(x_positions', 'x=%.1f')); {'Analytical'}], 'Location', 'eastoutside');

% % 2. Plot normalized profiles
% subplot(2,2,2);
% hold on;
% grid on;
% 
% for i = 1:length(x_positions)
%     x_pos = x_positions(i);
%     [~, idx] = min(abs(xc - x_pos));
%     u_normalized = u(idx,:) / max(u(idx,:));
%     plot(u_normalized, yc, 'Color', colors(i,:), 'LineWidth', 1.5);
% end
% 
% % Add normalized analytical
% % profile
% plot(u_analytical/max(u_analytical), y_analytical, 'k--', 'LineWidth', 2);
% 
% title('Normalized U-velocity Profiles');
% xlabel('u/u_{max}');
% ylabel('y-coordinate');

%% Create 3D visualization of velocity profiles
% figure(14);
% set(gcf, 'Position', [100 100 1200 700]);
% 
% % Number of profiles to show
% num_profiles = 20;
% x_indices = round(linspace(1, Nx, num_profiles));
% colors = jet(num_profiles);
% 
% % Create waterfall plot
% hold on;
% 
% for i = 1:num_profiles
%     idx = x_indices(i);
%     x_pos = xc(idx);
%     plot3(u(idx,:), ones(size(yc))*x_pos, yc, 'Color', colors(i,:), 'LineWidth', 2);
% end
% 
% % Add analytical profiles
% for i = 1:5:num_profiles
%     idx = x_indices(i);
%     x_pos = xc(idx);
% 
%     if x_pos >= Le_theoretical
%         u_analytical_at_x = 6 * u_in * y_analytical .* (H - y_analytical) / H^2;
%     else
%         blend_factor = x_pos / Le_theoretical;
%         u_fd = 6 * u_in * y_analytical .* (H - y_analytical) / H^2;
%         u_analytical_at_x = u_in * (1 - blend_factor) + u_fd * blend_factor;
%     end
% 
%     plot3(u_analytical_at_x, ones(size(y_analytical))*x_pos, y_analytical, 'k--', 'LineWidth', 2);
% end
% 
% view(30, 30);
% grid on;
% xlabel('u-velocity');
% ylabel('x-position');
% zlabel('y-position');
% title('3D Visualization of U-velocity Profile Development');

% Add colorbar to indicate x-position
% colormap(jet);
% cb = colorbar;
% cb.Ticks = linspace(0, 1, 5);
% cb.TickLabels = cellstr(num2str(linspace(0, L, 5)', '%.1f'));
% cb.Label.String = 'x-position';

% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%% Error analysis
%% Quantitative error analysis
% % Define the fully developed region
% fd_indices = find(xc >= entry_length);
% 
% % Calculate error metrics
% u_errors = zeros(length(fd_indices), Ny);
% l1_errors = zeros(length(fd_indices), 1);
% l2_errors = zeros(length(fd_indices), 1);
% linf_errors = zeros(length(fd_indices), 1);
% 
% for i = 1:length(fd_indices)
%     idx = fd_indices(i);
% 
%     % Calculate analytical solution
%     u_anal = 6 * u_in * yc .* (H - yc) / H^2;
% 
%     % Calculate errors
%     u_errors(i,:) = u(idx,:) - u_anal;
%     l1_errors(i) = mean(abs(u_errors(i,:)));
%     l2_errors(i) = sqrt(mean(u_errors(i,:).^2));
%     linf_errors(i) = max(abs(u_errors(i,:)));
% end
% 
% % Print error statistics
% fprintf('\n=== Error Analysis in Fully Developed Region ===\n');
% fprintf('Average L1 error: %.6e\n', mean(l1_errors));
% fprintf('Average L2 error: %.6e\n', mean(l2_errors));
% fprintf('Average L∞ error: %.6e\n', mean(linf_errors));
% 
% % Plot error norms
% figure(15);
% plot(xc(fd_indices), l1_errors, 'b-o', 'LineWidth', 1.5);
% hold on;
% plot(xc(fd_indices), l2_errors, 'r-s', 'LineWidth', 1.5);
% plot(xc(fd_indices), linf_errors, 'g-^', 'LineWidth', 1.5);
% grid on;
% title('Error Norms in Fully Developed Region');
% xlabel('x-coordinate');
% ylabel('Error Magnitude');
% legend({'L^1 Error', 'L^2 Error', 'L^∞ Error'}, 'Location', 'best');

%% 
%%%%%%%%%%%%%%%%% new structure mesh x= 10 units
% MATLAB Code for Structured Mesh Visualization
% For 2D Channel Flow (Plane Poiseuille flow)
% Channel dimensions: Height H = 1 unit, Length L = 10 units

% Problem Parameters
L = 10;              % Channel length [units]
H = 1;               % Channel height [units]
Nx = 100;            % Number of cells in x-direction
Ny = 20;             % Number of cells in y-direction

% Grid Generation
% Cell sizes
dx = L/Nx;
dy = H/Ny;

% Node coordinates (for mesh lines)
x = linspace(0, L, Nx+1);
y = linspace(0, H, Ny+1);

% Create figure for mesh visualization
figure;
hold on;

% Plot vertical grid lines (constant x)
for i = 1:length(x)
    plot([x(i) x(i)], [y(1) y(end)], 'k');
end

% Plot horizontal grid lines (constant y)
for j = 1:length(y)
    plot([x(1) x(end)], [y(j) y(j)], 'k');
end

% Set proper axis limits to show full domain length
axis([0 L 0 H]);

% Add labels and title
xlabel('x');
ylabel('y');
title('Structured Mesh');

% Optional: Add text to show domain dimensions
text(L/2, H*1.05, sprintf('L = %g units', L), 'HorizontalAlignment', 'center');
text(L*1.05, H/2, sprintf('H = %g unit', H), 'HorizontalAlignment', 'center', 'Rotation', 90);

% Turn off additional grid that might interfere with mesh display
grid off;
hold off;

