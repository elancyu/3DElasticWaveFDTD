clear all; clc;
% generate rough nanomesh input for 2D FDTD
% Define the maximum range
Lx = 600e-9;
Ly = 600e-9;
Lz = 145e-9;
Dp = 200e-9;
dx = 5.43e-9;
NX = round(Lx / dx);
NY = round(Lx / dx);
Nz = round(Lz / dx);
% define the pore size
nx = round(Dp / dx);
ny = round(Dp / dx);
hny = round(Dp * 0.5 / dx);
% define roughness level in the pore
nr = 3;
CX = (NX-1) / 2;
CY = (NY-1) / 2;
GEO = zeros(2*NX,NY,Nz+2);
for k = 2:Nz+1
    % 1st pore
    rho = ones(2*NX,NY);
    R = ones(nx, ny);
    % four sides
    Rxm = nr*rand(1,ny);
    Rxp = nr*rand(1,ny);
    Rym = nr*rand(nx,1);
    Ryp = nr*rand(nx,1);
    XM = zeros(nr,ny);
    XP = zeros(nr,ny);
    YM = zeros(nx, nr);
    YP = zeros(nx, nr);
    C = zeros(nr,nr);   % corner matrix
    for j = 1:ny
        nxm = Rxm(j);
        nxp = Rxp(j);
        for i = 1: nxm
            XM(i,j) = 1;
        end
        for i = 1: nxp
            XP(i,j) = 1;
        end
    end
    
    for i = 1: nx
        nym = Rym(i);
        nyp = Ryp(i);     
        for j = 1:nym
            YM(i,j) = 1;
        end
        for j = 1:nyp
            YP(i,j) = 1;
        end
    end
    
    XM = 1 - XM;
    YM = 1 - YM;
    YM = [C;YM;C];
    YP = [C;YP;C];
    R = [XM;R;XP];
    R = [YM R YP];
    ixstart = round(CX - nr - 0.5*nx);
    ixend = ixstart + 2*nr + nx - 1;
    iystart = round(CY - nr - 0.5*ny);
    iyend = iystart + 2*nr + ny - 1;
    rho(ixstart:ixend, iystart:iyend) = rho(ixstart:ixend, iystart:iyend)-R;
    % 2nd half
    R = ones(nx, hny);
    Rxm = nr*rand(1,hny);
    Rxp = nr*rand(1,hny);
    Ryp = nr*rand(nx,1);
    XM = zeros(nr,hny);
    XP = zeros(nr,hny);
    YP = zeros(nx, nr);
    C = zeros(nr,nr);   % corner matrix
    
    for j = 1:hny
        nxm = Rxm(j);
        nxp = Rxp(j);
        for i = 1: nxm
            XM(i,j) = 1;
        end
        for i = 1: nxp
            XP(i,j) = 1;
        end
    end
    
    for i = 1: nx
        nyp = Ryp(i);     
        for j = 1:nyp
            YP(i,j) = 1;
        end
    end
    
    XM = 1 - XM;
    YP = [C;YP;C];
    R = [XM;R;XP];
    R = [R YP];
    ixstart = NX + round(CX - nr - 0.5*nx);
    ixend = ixstart + 2*nr + nx - 1;
    iystart = 1;
    iyend = iystart + nr + hny - 1;
    rho(ixstart:ixend, iystart:iyend) = rho(ixstart:ixend, iystart:iyend)-R; 
    % 2nd other half
    R = ones(nx, hny);
    Rxm = nr*rand(1,hny);
    Rxp = nr*rand(1,hny);
    Rym = nr*rand(nx,1);
    XM = zeros(nr,hny);
    XP = zeros(nr,hny);
    YM = zeros(nx, nr);
    C = zeros(nr,nr);   % corner matrix
    
    for j = 1:hny
        nxm = Rxm(j);
        nxp = Rxp(j);
        for i = 1: nxm
            XM(i,j) = 1;
        end
        for i = 1: nxp
            XP(i,j) = 1;
        end
    end
    
    for i = 1: nx
        nym = Rym(i);     
        for j = 1:nym
            YM(i,j) = 1;
        end
    end
    
    XM = 1 - XM;
    YM = 1 - YM;
    YM = [C;YM;C];
    R = [XM;R;XP];
    R = [YM R];
    ixstart = NX + round(CX - nr - 0.5*nx);
    ixend = ixstart + 2*nr + nx - 1;
    iyend = NY;
    iystart = iyend + 1 - nr - hny;
    rho(ixstart:ixend, iystart:iyend) = rho(ixstart:ixend, iystart:iyend)-R;
    imagesc(rho);
    axis equal
    GEO(:,:,k) = rho;
end
fp = fopen('RNM.in','a+');
fprintf(fp,'%d %d %d\n', NX, NY, Nz+2);
fclose(fp);

for k = 1:Nz+2
    geo = GEO(:,:,k);
    save('RNM.in','geo','-ascii','-append');
end