% code to simulate uniform permeability model (with nonlinear, saturating
% permeability kinetics
function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = nonlinearPsims(options)
%% set up default simulation parameters
opt = struct();
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.vel = 1; % mitochondria velocity
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]

% in these simulations
% opt.kg = kg/(pi r^2 Km delta) where kg is glucose turnover per sec per
% mito

opt.kg = 0.1; % rate of glucose consumption (if linear) 
opt.Km = 1; % Michaelis-menten constant for glucose consumption
opt.PKm = 3; % Michaelis-menten constant for permeabililty

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 75; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-2; % time-step
opt.nstep = 1e5; % number of steps to run

%External glucose profile parameters
%default is linear external glucose profile
opt.c0 = 1;
opt.cend = 0.1; 
%Permeability term
opt.P = 0.1;
% tolerance for "small time derivative"
opt.dttol = 1e-4;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1;
opt.restart = 1; % flag to enable continuing previous sims

%%
% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end
% set up dimensionless parameters
%Nondimensionalize by L, L^2/D and Km
tscale = (opt.L)^2 / opt.D;
lscale = opt.L;
cscale = opt.Km;
Lh = opt.L/lscale;
kwh = opt.kw*tscale;
ksh = opt.ks*cscale*tscale; 
Dh = opt.D*tscale/lscale^2;
kgh = opt.kg*tscale;
Kmh = opt.Km/cscale;
PKmh = opt.PKm/cscale;
c0h = opt.c0/cscale;
cendh = opt.cend/cscale;
msizeh = opt.msize/lscale;
Ph = opt.P * tscale;
% spatial resolution
dx = Lh/(opt.gpts - 1);
delth = opt.delt/tscale;

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria 

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh*msizeh)); %lambda-hat
%% -------------
% initialize using analytical solution for constant linear consumption
kh = kgh*opt.nmito*msizeh/Lh;
a = (cendh-c0h)/Lh;
b = (cendh+c0h)/2;
h = Ph / (kh+Ph) * a;
j = Ph / (kh+Ph) * b;
Beta = sqrt((kh+Ph)/Dh);
if (Beta>100)
    gluc_calc_lin =  h*(xpos-Lh/2) + j;    
else
    A = h / (2*Beta*cosh(Beta*Lh/2));
    gluc_calc_lin = -2*A*sinh(Beta*(xpos-Lh/2)) + h*(xpos-Lh/2) + j;   
end

Beta = sqrt(Ph/Dh);
if (Beta>100)
    gluc_calc_max = a*(xpos-Lh/2)+b - kh*Kmh/Ph;
else
    A = a/(2*Beta*cosh(Beta*Lh/2));
    gluc_calc_max = -2*A*sinh(Beta*(xpos-Lh/2)) + a*(xpos-Lh/2)+b - kh*Kmh/Ph;
end

if (mean(gluc_calc_lin)>Kmh && min(gluc_calc_max)>0)
    % case with maxed out consumption
    gluc_calc = gluc_calc_max;
else
    gluc_calc=gluc_calc_lin;
end

%%
gluc_init = gluc_calc;
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
dtg = zeros(opt.gpts,1);
%define external glucose profile (nondimensionalized)
C_out = ((cendh - c0h) * (xpos - Lh/2) / Lh) + (cendh + c0h)/2;

ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 
normdtg = inf;
dtcutoff = opt.dttol*Kmh*kgh;
%tscale/cscale; %nondimensionalizing dttol
initglucint = dx * (sum(gluc_init) - (gluc_init(1)+gluc_init(end))/2);
%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    %Calculate distribution of total number of mitochondria
    ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
    ksx_int = dx * (sum(ksx)-(ksx(1)+ksx(end))/2);
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));
    Smito_int = dx * (sum(Smito)-(Smito(1)+Smito(end))/2);
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * Kmh * opt.nmito * msizeh) * (gluc(2:end-1) .* Tmito(2:end-1)) ./ (Kmh + gluc(2:end-1))...
        + ((Ph*PKmh*(C_out(2:end-1)-gluc(2:end-1)))./(PKmh +(C_out(2:end-1)-gluc(2:end-1))));
    normdtg = norm(dtg);
    gluc = gluc+dtg*delth;
    %implement reflecting boundary condition 
    gluc(1) = gluc(2);
    gluc(end) = gluc(end-1);
       
    if (any(gluc < -1e-3))
         disp('Concentration went negative. Try smaller timestep.')
        ftc = 1;
        return
    end
    step = step+1;
    
    if (step>opt.nstep)
        disp('Failed to converge')
        ftc = 1;
        return
    end
    
    if (opt.dodisplay && mod(step,opt.showevery)==0)
        
        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')
        %plot(xpos,C_out,'g--')
        plot(xpos,Tmito*initglucint,'r.-')
        title(sprintf('Step %d, normdtg: %f', step, normdtg))
        hold off
        drawnow
        
    end
end


        if (opt.dodisplay)
        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')
        plot(xpos,C_out,'g--')
        plot(xpos,Tmito*initglucint,'r.-')
        title(sprintf('Step %d, normdtg: %f', step, normdtg))
        hold off
        drawnow
        end
        
    
end
 