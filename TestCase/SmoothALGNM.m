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
NX = round( Lx / dx);
NY = round( Ly / dx);
NZ = round(Lz / dx) + 2;            % additional two layers for free boundary

nx = round(Dx / dx);
ny = round(Dy / dx);

% First Initialize with Zeros
RHO = zeros(NX, NY, NZ);

for k = 2:NZ - 1
    % Generate the geo layer by layer
    Map = ones(NX, NY);
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

fp = fopen('ALGNM.in','a+');
fprintf(fp,'%d %d %d\n', NX, NY, NZ);
fclose(fp);

for k = 1:NZ
    rho = RHO(:,:,k);
    save('ALGNM.in','rho','-ascii','-append');
end