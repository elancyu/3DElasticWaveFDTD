% This script generates the 3D smooth nanomesh
% assuming periodic boundary condition in X and Y directions
% simulation domain size
Lx = 34e-9;
Ly = 34e-9;
Lz = 22e-9;
Dp = 11e-9;
dx = 5.43e-10;          % grid spacing
Nx = round(Lx / dx);
Ny = round(Ly / dx);
Nz = round(Lz / dx);
Np = round(Dp / dx);

NM = zeros(Nx, Ny, Nz+2);

for i = 2:Nz+1
    NM2 = ones(Nx,Ny);
    ixs = round((Nx-Np)/2);
    ixe = round((Nx+Np)/2);
    iys = round((Ny-Np)/2);
    iye = round((Ny+Np)/2);
    nx = ixe - ixs;
    ny = iye - iys;
    PORE = zeros(nx, ny);
    NM2(ixs+1:ixe,iys+1:iye) = PORE;
    NM(:,:,i)=NM2;
end

for i = 1:Nz+2
    NMi = NM(:,:,i);
    save('SNM.in','NMi','-ascii','-append');
end