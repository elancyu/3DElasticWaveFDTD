% This script generates the 3D smooth nanomesh
% Size of the unit-cell
Lx = 600e-9;
Ly = 600e-9;
Lz = 145e-9;

% size of the pore.
Dx = 200e-9;
Dy = 200e-9;

% mesh grid size
dx = 5e-9;

% Numbers for defining the structure.
NX = round( 2* Lx / dx);
NY = round( Ly / dx);
NZ = round(Lz / dx) + 2;            % additional two layers for free boundary

nx = round(Dx / dx);
ny = round(Dy / dx);
hx = round(0.5 * Dx / dx);
hy = round(0.5 * Dy / dx);

% First Initialize with Zeros
RHO = zeros(NX, NY, NZ);

for k = 2:NZ - 1
    % Generate the geo layer by layer
    Map = ones(NX, NY);
    % Voids at the corners
    Map(1:hx, 1:hy) = zeros(hx, hy);
    Map(1:hx, NY - hy + 1: NY) = zeros(hx, hy);
    Map(NX - hx + 1:NX, 1:hy) = zeros(hx, hy);
    Map(NX - hx + 1: NX, NY - hy + 1: NY) = zeros(hx, hy);
    % The central void
    startx = round(0.5 * (NX - nx)) + 1;
    endx = startx + nx - 1;
    starty = round(0.5 * (NY - ny)) + 1;
    endy = starty + ny - 1;
    Map(startx: endx, starty: endy) = zeros(nx, ny);
    %imagesc(Map);
    %axis equal;
    RHO(:,:,k)= Map;
end

fp = fopen('STGNM.in','a+');
fprintf(fp,'%d %d %d\n', NX, NY, NZ);
fclose(fp);
for k = 1:NZ
    rho = RHO(:,:,k);
    save('STGNM.in','rho','-ascii','-append');
end