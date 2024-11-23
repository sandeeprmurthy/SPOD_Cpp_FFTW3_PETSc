clear all;
clc;
close all;

addpath(genpath('./lib'));
addpath(genpath('./lib2'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the data path and file name
disp('Read the file jetLES.h5 ...')
data_path = '../'; % Replace with the actual path
file_name = fullfile(data_path, 'jetLES.h5');

% Open the HDF5 file
data = h5read(file_name, '/data');    % Flow fields
grid = h5read(file_name, '/grid');    % Grid points
dt   = h5read(file_name, '/dt');      % Unit in seconds
dt   = dt(1);                         % Extract the first value if needed

% Calculate dimensions
ng   = size(grid, 1);                 % Number of grid points
nt   = size(data, 1);                 % Number of snapshots
nx   = size(data, 2);                 % Total grid points * number of variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Read the vector q ...')
q = PetscBinaryRead('./spod_modes/spod_mode_freq_5_mode_0.dat', 'indices','int32','precision','float64','complex',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qlevels = -0.06:0.006:0.06;    % Contour levels
qname = '$\phi_p$';            % Label for the colorbar

% Call jet_contour function
disp('Plot the eigenvector ...')
% Create figure
fig = figure('Units', 'inches', 'Position', [0, 0, 12, 8]);

% Perform Delaunay triangulation
x = grid(1,:);
y = grid(2,:);
tri = delaunay(x, y);

% Plot filled triangular contours
trisurf(tri, x, y, real(q), 'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal
view(2); % Set view to 2D
shading interp; % Smooth shading
colorbar;
xlabel('X-axis');
ylabel('Y-axis');    

% Save the figure if needed
save_fig = true;
if save_fig
    save_path = './'; % Replace with your save path
    saveas(fig, fullfile(save_path, 'SPOD_mode.png'));
end

