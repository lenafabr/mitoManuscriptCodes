function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options)
%% set up default simulation parameters

opt = struct();

opt.kg = 1; % rate of glucose consumption
opt.c0 = 0.01; % fixed glucose concentration

opt.L = 250; % domain size
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]
opt.Km = 0.1; % Michaelis constant for glucose consumption

opt.avgmitoden=0.3; % mitochondrial density, #/length^3

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

% use the infinite ks limit
opt.highks = 0;

opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-5; % time-step
opt.nstep = 1e6; % number of steps to run

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% tolerance for "small time derivative"
opt.dttol = 1e-3;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1000;

opt.restart = 1; % flag to enable continuing previous sims

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

% set up dimensionless parameters
%Nondimensionalize by L, L^2/D and L^3
tscale = (opt.L)^2 / opt.D;
lscale = opt.L;

Lh = opt.L/lscale;
kwh = opt.kw*tscale;
ksh = opt.ks*tscale;
Dh = opt.D*tscale/lscale^2;
kgh = opt.kg*tscale;
c0h = opt.c0*lscale^3;
Kmh = opt.Km*lscale^3;
avgmitodenh = opt.avgmitoden*lscale^3;

% spatial resolution
dx = Lh/(opt.gpts - 1);


%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh*Kmh/(kgh*avgmitodenh*Lh^2)); %lambda-hat
if (~isempty(opt.startgluc))
    gluc_init=opt.startgluc;
else
    Vmax = kgh*avgmitodenh;
    gluc_init = Vmax/(2*Dh)*(xpos-Lh/2).^2 + c0h - Vmax*Lh^2/(8*Dh);
    if (min(gluc_init)<0)
        gluc_init = c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
    end    
    %if (
end
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
dtg = zeros(opt.gpts,1);

ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 

normdtg = inf;

dtcutoff = opt.dttol;

%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    
    %Calculate distribution of total number of mitochondria
    ksx = ksh * gluc ./ (Kmh + gluc);    
    ksx_int = dx/2 * (ksx(1) + ksx(end) + 2*sum(ksx(2:end-1)));
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
     if (opt.highks)
         Tmito = (ksx/kwh) ./ (ksx_int/kwh);
     end
    
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));    
    Smito_int = dx/2 * (Smito(1) + Smito(end) + 2*sum(Smito(2:end-1)));
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * avgmitodenh*Lh) * (gluc(2:end-1) .* Tmito(2:end-1)) ./ (Kmh + gluc(2:end-1));
    normdtg = norm(dtg)/(kgh * avgmitodenh);

    gluc = gluc+dtg*opt.delt;
       
    if (any(isnan(gluc)))
        error('problem with gluc going to NaN')
    end
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
        % display distributions as you go
        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')        
        plot(xpos,Tmito*c0h*Lh,'r.-')
        hold off
        legend('orig','gluc','mito')
        title(sprintf('Norm dg/dt: %f',normdtg))
        drawnow
    end
end

end

