classdef PML < handle

    properties
        %%% Domain info
        xmin, ymin % coordinate of bottom left corner (real number)
        Lx, Ly % domain physical size (Euler + PML) (real number)
        Nx, Ny % Grid resolution (integer)
        Ix, Iy, I

        %%% PML info
        Px, Py, % width of PML in number of grids (integer)

        sigmax % damping coeff
        sigmay % damping coeff
        sxarr % damping coeff array
        syarr % damping coeff array
        power % damping profile power

        %%% Calculated variables
        dx, dy % grid size (real number)
        dt % time step size (real number)
        X, Y % grid coordinate (IND -> R)
        NN % total number of grid points (2 * PML and domain for each dimension)


        %%% Other info

        % Field and physical variables
        rho
        u
        v     
        p
        q3
        q4
        Mach
        k

        % Storage
        rhostorage
        ustorage
        vstorage
        pstorage
        q3storage
        q4storage
        krhostorage
        kustorage
        kvstorage
        kpstorage

        % DRP and RK coefficients
        b0, b1, b2, b3
%         a0, a1, a2, a3
        a

        % Misc
        xlen, ylen
        timelevel
        tempvar
        tempvar2

   

    end
    methods
        %% Constructor
        function pml = PML(varargin)
            switch nargin
                case 0
                    return
                case 8
                    pml.xmin = varargin{1};
                    pml.ymin = varargin{2};
                    pml.Lx = varargin{3};
                    pml.Ly = varargin{4};
                    pml.Nx = varargin{5};
                    pml.Ny = varargin{6};
                    pml.Px = varargin{7};
                    pml.Py = varargin{8};
                    pml.build();
                otherwise
                    error('PML constructor inputs: xmin,ymin,Lx,Ly,Nx,Ny,Px,Py')
            end
        end
        %%% Build/Initialiser
        function build(pml)
            % build - helper function for the constructor
            % given inputs, computer intialise all other class properties
            pml.dx = pml.Lx/pml.Nx;
            pml.dy = pml.Ly/pml.Ny;
            pml.xlen = pml.Nx; % Total grid points of x dimension
            pml.ylen = pml.Ny; % "" for y
            pml.dt=(min([pml.dx,pml.dy])/2)/10;
            [pml.Ix,pml.Iy] = ndgrid(1:pml.Nx,1:pml.Ny); % Linear to sub indexes
            pml.Ix = pml.Ix(:);
            pml.Iy = pml.Iy(:);
            pml.X = pml.xmin + (pml.Ix-1)*pml.dx;
            pml.Y = pml.ymin + (pml.Iy-1)*pml.dy;
            pml.I  = pml.sub2ind(pml.Ix,pml.Iy);
            pml.NN = pml.Nx * pml.Ny;
%             pml.setupPML()



            
            % Optimised time marching coefficients
            pml.b0 = 2.3025580888383;
            pml.b1 = -2.4910075998482;
            pml.b2 = 1.5743409331815;
            pml.b3 = -0.3858914221716;
            % Optimised DRP coefficients
%             pml.a0 = 0;
%             pml.a1 = 0.7708823805182255; % Also a1=a-1, etc
%             pml.a2 = - 0.1667059044145804;
%             pml.a3 = 0.0208431427703117;
            % I hate this so much as a(1) = a0, but whatever FML
            pml.a(1) = 0; % a0
            pml.a(2) = 0.7708823805182255; % a1, Also a1=a-1, etc
            pml.a(3) = - 0.1667059044145804; % a2
            pml.a(4) = 0.0208431427703117; % a3

            % Storage for DRP and multilevel time marching
            % First row = n, second row = n-1, third row = n-2, fourth
            % row = n-3
            pml.rhostorage = zeros(4,pml.NN);
            pml.ustorage = pml.rhostorage;
            pml.vstorage = pml.rhostorage;
            pml.pstorage = pml.rhostorage;
            pml.q3storage = pml.rhostorage;
            pml.q4storage = pml.rhostorage;
            pml.krhostorage = zeros(4,pml.NN);
            pml.kustorage = pml.krhostorage;
            pml.kvstorage = pml.krhostorage;
            pml.kpstorage = pml.krhostorage;


        end


        function setupPML(pml)
            % setupPML - creates Perfectly Matched Layers within the
            % domain, with profiles determined by a power law with a max
            % sigmax and sigmay
            % Requires initialisation in setup script of: sigmax, sigmay,
            % power

            pml.tempvar = pml.sigmax * ( (((1:pml.Px)*pml.dx)/(pml.Px * pml.dx)) .^ pml.power); % Power law. Damping coefficient profile in PML layer
            pml.sxarr = zeros(size(pml.I)); % Initialise array for entire domain (x dimension)
            pml.sxarr(pml.Ix < (pml.Px+1)) = repmat(pml.tempvar(end:-1:1),1,pml.Ny)'; % Place PML damping coefficient values in domain (start of x dimension i.e. LHS)
            pml.sxarr(pml.Ix > pml.Nx-pml.Px) = repmat(pml.tempvar,1,pml.Ny)'; % Place PML damping coefficient values in domain (end of x dimension i.e. RHS)
            pml.sxarr = pml.sxarr';

            pml.tempvar = pml.sigmay * ( (((1:pml.Py)*pml.dy)/(pml.Py * pml.dy)) .^ pml.power); % Power law
            pml.syarr = zeros(size(pml.I));
            pml.tempvar2 = repmat(pml.tempvar, pml.Nx,1);
            pml.tempvar2 = pml.tempvar2(:);
            pml.syarr(pml.Iy < (pml.Py+1)) = pml.tempvar2(end:-1:1);
            pml.syarr(pml.Iy > pml.Ny-pml.Py) = pml.tempvar2;


        end

        %%% Derivative functions

        function [deriv] = drhodx(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((l + j) < 1) || ((l + j) > pml.xlen) % Ghost node catch
                    j = 0; % If ghost node, set finite diff index to 0 (ghost node equal to inner mesh value)
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.rhostorage(1,pml.xlen*(m-1)+l+j);
            end
            deriv = (1/pml.dx) * deriv;
        end
        function [deriv] = dudx(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((l + j) < 1) || ((l + j) > pml.xlen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.ustorage(1,pml.xlen*(m-1)+l+j);
            end
            deriv = (1/pml.dx) * deriv;
        end
        function [deriv] = dvdx(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((l + j) < 1) || ((l + j) > pml.xlen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.vstorage(1,pml.xlen*(m-1)+l+j);
            end
            deriv = (1/pml.dx) * deriv;
        end
        function [deriv] = dpdx(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((l + j) < 1) || ((l + j) > pml.xlen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.pstorage(1,pml.xlen*(m-1)+l+j);
            end
            deriv = (1/pml.dx) * deriv;
        end
        function [deriv] = dpdy(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((m + j) < 1) || ((m + j) > pml.ylen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.pstorage(1,pml.xlen*(m+j-1)+l);
            end
            deriv = (1/pml.dy) * deriv;
        end
        function [deriv] = dvdy(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((m + j) < 1) || ((m + j) > pml.ylen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.vstorage(1,pml.xlen*(m+j-1)+l);
            end
            deriv = (1/pml.dy) * deriv;
        end
        function [deriv] = dq3dy(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((m + j) < 1) || ((m + j) > pml.ylen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.q3storage(1,pml.xlen*(m+j-1)+l);
            end
            deriv = (1/pml.dy) * deriv;
        end
        function [deriv] = dq4dy(pml, l, m)
            deriv = 0;
            for j=-3:3
                if ((m + j) < 1) || ((m + j) > pml.ylen)
                    j = 0;
                end
                deriv = deriv +  pml.a(abs(j)+1) * pml.q4storage(1,pml.xlen*(m+j-1)+l);
            end
            deriv = (1/pml.dy) * deriv;
        end



        function DRPStep(pml, timelevel)
        % DRPStep - updates field variables by one step of DRP and RK
        % Solves ODE for one time step and updates values



        % Check if first iteration, if so assign initial values to storage
        % variables
        if timelevel == 1
            pml.rhostorage(1,:) = pml.rho;
            pml.ustorage(1,:) = pml.u;
            pml.vstorage(1,:) = pml.v;
            pml.pstorage(1,:) = pml.p;
            pml.q3storage(1,:) = 0;
            pml.q4storage(1,:) = 0;
        end

        % Main loop. At each grid point, calculate intermediate k and then time derivative to get field variables at time n+1. Ghost nodes handled in DRP
        % functions
        for m=1:pml.ylen
            for l=1:pml.xlen
                idx = pml.xlen*(m-1)+l;
                pml.krhostorage(1,idx) = - pml.Mach*pml.drhodx(l,m)...
                    - pml.dudx(l,m) - pml.dvdy(l,m)...
                    - pml.sxarr(idx)*pml.dq3dy(l,m)...
                    - pml.sxarr(idx)*pml.rho(idx)...
                    - (pml.sxarr(idx)*pml.Mach / (1-pml.Mach^2))...
                    * (pml.Mach*pml.rho(idx) + pml.u(idx)); % Spatial discretisation for rho
                pml.rho(idx) = pml.rhostorage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.krhostorage(1,idx)...
                    + pml.b1 * pml.krhostorage(2,idx)...
                    + pml.b2 * pml.krhostorage(3,idx)...
                    + pml.b3 * pml.krhostorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)
                pml.kustorage(1,idx) = - pml.Mach*pml.dudx(l,m)...
                    - pml.dpdx(l,m) - pml.sxarr(idx)*pml.u(idx)...
                    - (pml.sxarr(idx)*pml.Mach / (1-pml.Mach^2))...
                    * (pml.Mach*pml.u(idx) + pml.p(idx)); % Spatial discretisation for u
                pml.u(idx) = pml.ustorage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.kustorage(1,idx)...
                    + pml.b1 * pml.kustorage(2,idx)...
                    + pml.b2 * pml.kustorage(3,idx)...
                    + pml.b3 * pml.kustorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)
                pml.kvstorage(1,idx) = - pml.Mach*pml.dvdx(l,m)...
                    - pml.dpdy(l,m) - pml.sxarr(idx)*pml.dq4dy(l,m)...
                    - pml.sxarr(idx)*pml.v(idx)...
                    - (pml.sxarr(idx)*pml.Mach^2 / (1-pml.Mach^2))*pml.v(idx); % Spatial discretisation for v
                pml.v(idx) = pml.vstorage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.kvstorage(1,idx)...
                    + pml.b1 * pml.kvstorage(2,idx)...
                    + pml.b2 * pml.kvstorage(3,idx)...
                    + pml.b3 * pml.kvstorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)
                pml.kpstorage(1,idx) = - pml.Mach*pml.dpdx(l,m)...
                    - pml.dudx(l,m) - pml.dvdy(l,m)...
                    - pml.sxarr(idx)*pml.dq3dy(l,m)...
                    - pml.sxarr(idx)*pml.rho(idx)...
                    - (pml.sxarr(idx)*pml.Mach / (1-pml.Mach^2))...
                    * (pml.Mach*pml.p(idx) + pml.u(idx)); % Spatial discretisation for p
                pml.p(idx) = pml.pstorage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.kpstorage(1,idx)...
                    + pml.b1 * pml.kpstorage(2,idx)...
                    + pml.b2 * pml.kpstorage(3,idx)...
                    + pml.b3 * pml.kpstorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)
                pml.q3(idx) = pml.q3storage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.vstorage(1,idx)...
                    + pml.b1 * pml.vstorage(2,idx)...
                    + pml.b2 * pml.vstorage(3,idx)...
                    + pml.b3 * pml.vstorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)
                pml.q4(idx) = pml.q4storage(1,idx) + pml.dt...
                    * ( pml.b0 * pml.pstorage(1,idx)...
                    + pml.b1 * pml.pstorage(2,idx)...
                    + pml.b2 * pml.pstorage(3,idx)...
                    + pml.b3 * pml.pstorage(4,idx) ); % RK (simple 4step for now, will do LDDRK56 eventually)


            end

        end

        % Reassign storage for time marching scheme (shift down so row 1 is
        % assigned to row 2 etc. Then put current variable in at row 1 (for
        % time n))
        % 2N storage scheme

        pml.rhostorage = [pml.rho'; pml.rhostorage(1:end-1,:)];
        pml.ustorage = [pml.u'; pml.ustorage(1:end-1,:)];
        pml.vstorage = [pml.v'; pml.vstorage(1:end-1,:)];
        pml.pstorage = [pml.p'; pml.pstorage(1:end-1,:)];
        pml.q3storage = [pml.q3; pml.q3storage(1:end-1,:)];
        pml.q4storage = [pml.q4; pml.q4storage(1:end-1,:)];
        pml.krhostorage = [zeros(1,length(pml.krhostorage)); pml.krhostorage(1:end-1,:)];
        pml.kustorage = [zeros(1,length(pml.kustorage)); pml.kustorage(1:end-1,:)];
        pml.kvstorage = [zeros(1,length(pml.kvstorage)); pml.kvstorage(1:end-1,:)];
        pml.kpstorage = [zeros(1,length(pml.kpstorage)); pml.kpstorage(1:end-1,:)];

        end


        %%% Helper functions

        function ind = sub2ind(pml,subx,suby)
            % sub2ind - return index given the x y positions

            periodify = @(ind,N) mod(ind-1,N)+1;
            ind = sub2ind([pml.Nx,pml.Ny],...
                periodify(subx,pml.Nx),periodify(suby,pml.Ny));
        end


        function vv = vec2grid(pml,v)
            % converts linear array into 2D array

            vv = reshape(v, pml.Nx, pml.Ny);

        end
end
end