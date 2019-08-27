clear all; clc;
fprintf('Creating Staggered Nanomesh Geometry...\n');
% This script generates the 3D smooth nanomesh
% assuming periodic boundary condition in X and Y directions
% simulation domain size
Lx = 600e-9;
Ly = 600e-9;
Lz = 145e-9;
Dp = 200e-9;
dx = 5e-9;          % grid spacing
Nx = round(Lx / dx);
Ny = round(Ly / dx);
Nz = round(Lz / dx);
Np = round(Dp / dx);

NM = zeros(2 * Nx, Ny, Nz+2);

for i = 2:2
    NM2 = ones(2*Nx,Ny);
    ixs = round((Nx-Np)/2);
    ixe = round((Nx+Np)/2);
    iys = round((Ny-Np)/2);
    iye = round((Ny+Np)/2);
    nx = ixe - ixs;
    ny = iye - iys;
    PORE = zeros(nx, ny);
    NM2(ixs+1:ixe,iys+1:iye) = PORE;
    % 1st half pores
    ixs = Nx + round((Nx-Np)/2);
    ixe = Nx + round((Nx+Np)/2);
    iys = 0;
    iye = round(Np/2);
    nx = ixe - ixs;
    ny = iye - iys;
    PORE = zeros(nx, ny);
    NM2(ixs+1:ixe,iys+1:iye) = PORE;
    % 2nd half pores
    ixs = Nx + round((Nx-Np)/2);
    ixe = Nx + round((Nx+Np)/2);
    iys = round(Ny - Np/2);;
    iye = Ny;
    nx = ixe - ixs;
    ny = iye - iys;
    PORE = zeros(nx, ny);
    NM2(ixs+1:ixe,iys+1:iye) = PORE;
    % end
    NM(:,:,i)=NM2;
    imagesc(NM2);
    axis equal
    axis off
end

% output the size of the simulation domain
fp = fopen('SNM.in','a+');
fprintf(fp,'%d  %d  %d\n', Nx, Ny, Nz+2);
fclose(fp);

% output the matrices
for i = 1:Nz+2
    NMi = NM(:,:,i);
    save('SNM.in','NMi','-ascii','-append');
end