
[nx, ny, nz, matrix] = loadmat3d("build/cart.txt");
[Y, X, Z] = meshgrid(1:ny, 1:nx, 1:nz);



figure;
h = gcf;
slice(Y, X, Z, matrix, [ny/2 ny], [nx/2 nx], [1 nz/2]);
colormap(jet);
colorbar;
caxis([0 65]);
