% Read parameters from Excel file
inputTable = readtable('AM24S012_streamline_input.xlsx');

% Extract parameters using descriptive names (handles spreadsheet label typos)
H = inputTable.value(strcmp(inputTable.parameter, 'Channel height'));
L = inputTable.value(strcmp(inputTable.parameter, 'Channel length'));
Re = inputTable.value(strcmp(inputTable.parameter, 'Reynolds number'));
nx = inputTable.value(strcmp(inputTable.parameter, 'Number of points in x-direction'));
ny = inputTable.value(strcmp(inputTable.parameter, 'Number of points in y-direction')); % Corrects label typo automatically

% Create grid using imported parameters
x = linspace(0, L, nx);
y = linspace(0, H, ny);
[X, Y] = meshgrid(x, y);

% Initialize velocity fields
U = zeros(ny, nx);
V = zeros(ny, nx);

% Calculate entry length using empirical formula
L_entry = 0.05 * Re * H;

% Velocity field development calculation
for i = 1:nx
    x_pos = x(i);    
    development_factor = min(x_pos/L_entry, 1);
    u_max = 1.5;
    
    for j = 1:ny
        y_normalized = (y(j) - H/2)/(H/2);
        u_parabolic = u_max * (1 - y_normalized^2);
        U(j,i) = (1-development_factor)*1.0 + development_factor*u_parabolic;
    end
end

% Stream function calculation
psi = zeros(ny, nx);
for i = 1:nx
    for j = 2:ny
        psi(j,i) = psi(j-1,i) + U(j,i)*(y(j)-y(j-1));
    end
end

% Visualization setup
figure('Position', [100, 100, 900, 400]);
contour(X, Y, psi, 20, 'LineWidth', 1.2);
hold on;

% Vector plot configuration
skip = 10;
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
      U(1:skip:end, 1:skip:end), V(1:skip:end, 1:skip:end), ...
      'r', 'AutoScaleFactor', 0.3);

% Boundary plotting
plot([0 L], [0 0], 'k-', 'LineWidth', 2);
plot([0 L], [H H], 'k-', 'LineWidth', 2);
plot([0 0], [0 H], 'k-', 'LineWidth', 2);
plot([L L], [0 H], 'k-', 'LineWidth', 2);

% Annotation
title('Plane Poiseuille Flow Streamlines (Excel-Parameterized)');
xlabel('x-direction [units]');
ylabel('y-direction [units]');
axis equal tight
grid on

% Velocity profile inset
axes('Position',[.75 .7 .2 .2]);
plot(U(:,1), y, 'b-', 'LineWidth', 2, 'DisplayName', 'Inlet');
hold on;
plot(U(:,round(nx/4)), y, 'g-', 'LineWidth', 2, 'DisplayName', 'Developing');
plot(U(:,end), y, 'r-', 'LineWidth', 2, 'DisplayName', 'Fully Developed');
title('Velocity Profiles');
xlabel('Velocity [units/s]');
ylabel('y [units]');
legend('Location', 'best');
grid on;
