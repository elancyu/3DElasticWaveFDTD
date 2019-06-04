% This script generate the 3D smooth NW.
clear all;clc;
% define the geometry
dx = 5.43e-10;          % grid spacing
Lx = 10e-9;                  % Length in X direction.
Ly = 10e-9;                  % Length in Y direction.
Lz = 10e-9;                  % Length in Z direction: thickness

% Turn to integer
Nx = round(Lx / dx);
Ny = round(Ly / dx);
Nz = round(Lz / dx);

% NW Solid
NW = ones(Nx, Ny, Nz);

% Side Wall Padding
L = zeros(1, Ny, Nz);
R = zeros(1, Ny, Nz);

NW = [L; NW; R];

% Top/Bottom Padding
T = zeros(Nx, Ny, 1);
B = zeros(Nx, Ny, 1);

GEO = zeros(Nx+2, Ny, Nz + 2);
GEO(:,:,2:Nz+1) = NW;

% CHECK each layer
% for i = 1:Nz+2
%     figure(i)
%     imagesc(GEO(:,:,i));colormap jet;
%     title([num2str(i) ' layer']);
%     legend('X','Y');
% end

% output the geo layer by layer
for i = 1:Nz+2
    X = GEO(:,:,i);
    save('SNW.in','X','-ascii','-append');
end