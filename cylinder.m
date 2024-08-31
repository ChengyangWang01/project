[nr, ntheta, nz, T] = loadmat3d("build/cylinder.txt");

theta = linspace(0, 2*pi, ntheta); 
r     = linspace(0, 1.0/3.0, nr);
z     = linspace(0, 1, nz);

[R, Theta, Z] = ndgrid(r, theta, z);

X = R .* cos(Theta);
Y = R .* sin(Theta);

% figure;
% for k = 1:length(z)
%   sliceZ = Z(:,:,k);
%   sliceT = T(:,:,k);
%   sliceX = X(:,:,k);
%   sliceY = Y(:,:,k);

%   surf(sliceX, sliceY, sliceZ, sliceT, 'EdgeColor', 'none');
%   hold on;
% end


% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('3D Heat Conduction in Cylindrical Coordinates');
% colorbar;
% view(3);

% colormap("jet")
% axis tight;

theta_min = -pi / 4;
theta_max = 3*pi/2;


theta_idx = find(theta >= theta_min & theta <= theta_max);

% Filter theta and corresponding X, Y, T data in the angle range
theta_filtered = theta(theta_idx);
R_filtered = R(:, theta_idx, :);
Theta_filtered = Theta(:, theta_idx, :);
Z_filtered = Z(:, theta_idx, :);
X_filtered = X(:, theta_idx, :);
Y_filtered = Y(:, theta_idx, :);
T_filtered = T(:, theta_idx, :);

% plot figure
figure;
for k = 1:length(z)
  sliceZ = Z_filtered(:,:,k);
  sliceT = T_filtered(:,:,k);
  sliceX = X_filtered(:,:,k);
  sliceY = Y_filtered(:,:,k);

  surf(sliceX, sliceY, sliceZ, sliceT, 'EdgeColor', 'none');
  hold on;
end


xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Heat Conduction in Cylindrical Coordinates (Filtered by Angle)');
colorbar;
view(3);

colormap("jet")
axis tight;

% figure;
% h = gcf;
% slice(Y, X, Z, matrix, [ny/2 ny], [nx/2 nx], [1 nz/2]);
% colormap(jet);
% colorbar;
% caxis([0 65]);