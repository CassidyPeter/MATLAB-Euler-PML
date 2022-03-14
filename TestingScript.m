
rng(0)
clc; clear; close all;

% Domain constants and setup
xmin=0;
ymin=0;
Lx=3;
Ly=3;
Nx=100;
Ny=100;
Px=15;
Py=15;
T=0.5; % Sim time in seconds
% Mach = 0.5;
Mach = 0.2;

% Initialise PML object
pml = PML(xmin,ymin,Lx,Ly,Nx,Ny,Px,Py);

% Setup PML
pml.Mach = Mach;
pml.power = 2;
pml.sigmax = 1;
pml.sigmay = 1;
pml.setupPML();

% % Initial conditions - Gaussian pulse
% epsilon = 1; % perturbation amplitude
% r0 = 9 / log(2); % Characteristic dimension of pulse (variance)
% r0 = 0.1;
% x0 = 1.5; % Source position
% y0 = 1.5; % Source position
% pml.p = epsilon * exp(-(((pml.X-x0).^2 + (pml.Y-y0).^2)/r0));
% pml.p(pml.sxarr'~=0 | pml.syarr~=0) = 0;
% pml.rho = pml.p;
% pml.u = zeros(size(pml.I));
% pml.v = pml.u;


% Initial conditions - acoustic, vorticity, entropy pulses
epsilon = 1; % perturbation amplitude
r0 = 9 / log(2); % Characteristic dimension of pulse (variance)
r0 = 0.1;
x0 = 1.5; % Source position
y0 = 1.5; % Source position
pml.p = epsilon * exp(-(((pml.X-x0).^2 + (pml.Y-y0).^2)/r0));
pml.p(pml.sxarr'~=0 | pml.syarr~=0) = 0;
pml.rho = epsilon * exp(-(((pml.X-1).^2 + (pml.Y-y0).^2)/r0)) + epsilon * exp(-(((pml.X-2).^2 + (pml.Y-y0).^2)/r0));
pml.rho(pml.sxarr'~=0 | pml.syarr~=0) = 0;
pml.u = epsilon * pml.Y .* exp(-(((pml.X-x0).^2 + (pml.Y-y0).^2)/r0));
pml.u(pml.sxarr'~=0 | pml.syarr~=0) = 0;
pml.v = -epsilon * (pml.X-x0) .* exp(-(((pml.X-x0).^2 + (pml.Y-y0).^2)/r0));
pml.v(pml.sxarr'~=0 | pml.syarr~=0) = 0;



% Initialise video
myVideo = VideoWriter('PMLTestVideo'); % Video animation line
myVideo.FrameRate = 10; % Video animation line
open(myVideo) % Video animation line

% Initialise graphics
hfsol = figure;
hsol = surf(pml.vec2grid(pml.X),...
             pml.vec2grid(pml.Y),...
            double( pml.vec2grid(pml.u) ));
set(hsol,'edgecolor','none')
% [~,hsol] = contour(pml.vec2grid(pml.X),...
%              pml.vec2grid(pml.Y),...
%             double( pml.vec2grid(pml.rho) ));
cameratoolbar
axis([pml.xmin,pml.xmin+pml.Lx,pml.ymin,pml.ymin+pml.Ly,-1,1])
light
colorbar
drawnow

itermax = ceil(T/pml.dt);
tt = (1:itermax)*pml.dt;
errrec = nan(itermax,1);

for iter = 1:itermax
    pml.DRPStep(iter)
    set(hsol,'ZData',double(pml.vec2grid(pml.u)))
    drawnow
    pause(0.01) % Video animation line
    fprintf(['Iteration ' num2str(iter) ' out of ' num2str(itermax) '\n'])

    frame = getframe(gcf); % Video animation line
    writeVideo(myVideo,frame); % Video animation line

end

close(myVideo) % Video animation line