%Numerical simulations for when diffusive timescale is slow
function [gluc,M,gluc_init,opt,xpos,ftc,normdg,glucplus,glucminus] = slowdiffusionsim(options)
%% define options structure
opt.Km = 0.03e-3*6e23/1000/1e12; % Michaelis-menten constant for glucose consumption
opt.PKm = 2.87e-3*6e23/1000/1e12; % Michaelis-menten constant for permeabililty
opt.r = 0.4; % radius of axon
opt.ks = 19/opt.Km; % rate of stopping is ks*[gluc], high ks regime
opt.kg = 1.3e5; % rate of glucose consumption
opt.L = 1e3;
opt.D = 140;% glucose diffusion coefficient
opt.nmito = 0.3*opt.L*pi*opt.r^2; % number of mitochondria
opt.kw = 1; % rate of starting a walk
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-6; % time step
opt.nstep = 1e10; % number of steps to run
%default is linear external glucose profile
opt.c0 = 1e-3*6e23/1000/1e12; % Using upper limit on brain glucose
opt.cend = 0.1e-3*6e23/1000/1e12; % lower limit brain glucose
opt.avgmitoden=0.3; % mitochondrial density, #/length^3

opt.P = 2e-2; % permeability in microns/sec 
opt.Pbykg = opt.P/opt.kg; % ratio, will be provided during run loops
% tolerance for glucose difference
opt.dttol = 1; %gluc is expressed in terms of #/vol
% displaying plots
opt.dodisplay = 1;
opt.showevery = 1000;
opt.restart = 1; % flag to enable continuing previous sims

% high ks limit
opt.highks = 0;
%%
% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end


% set up dimensionless parameters
%Nondimensionalize by L, L^2/D 
tscale = (opt.L)^2 / opt.D; % time scale
lscale = opt.L; % length scale
cscale = (opt.L)^(-3); % concentration scale
Lh = opt.L/lscale;
kwh = opt.kw*tscale;
ksh = opt.ks*cscale*tscale; 
Kmh = opt.Km/cscale;
PKmh = opt.PKm/cscale;
c0h = opt.c0/cscale;
cendh = opt.cend/cscale;
Pbykgh = opt.Pbykg/lscale;
rh = opt.r/lscale;
% spatial resolution
dx = Lh/(opt.gpts - 1);
avgmitodenh = opt.avgmitoden*lscale^3;

xpos = linspace(0,Lh,opt.gpts)';
%% 
% iterate over gluc and mitodensity
C_out = ((cendh - c0h) * (xpos - Lh/2) / Lh) + (cendh + c0h)/2;
gluc_init = ((cendh - c0h) * (xpos - Lh/2) / Lh) + (cendh + c0h)/2;
gluc = gluc_init;
ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
ksx_int = dx * (sum(ksx)-(ksx(1)+ksx(end))/2);
% if (opt.highks)
%     M  = avgmitodenh  * (ksx/kwh) ./ ((ksx_int/kwh));
% else
%     M  = avgmitodenh  * (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
% end
M  = avgmitodenh * (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
step = 0;
ftc = 0;
normdg = inf;
while (normdg/mean(gluc) > opt.dttol)
    glucold = gluc;
    alpha = 2 * Pbykgh * PKmh ./ (rh * M);
    a = 1 - alpha;
    b = ((alpha-1) .* C_out) - alpha*Kmh -PKmh;
    c = alpha * Kmh .* C_out;
    glucplus = (-b + sqrt(b.^2 - 4*a.*c)) ./ (2*a); %adding positive sqrt
    glucminus = (-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a); %adding negative root
 
% Assign the value of gluc which is positive everywhere
    plusgood = all(glucplus>-1e-3 & glucplus < C_out+1e-3);
    minusgood = all(glucminus>-1e-3 & glucminus < C_out+1e-3);

    if (plusgood && minusgood)
        error('two possible solutions')
    elseif (plusgood)
        gluc = (glucplus+gluc)/2;
    elseif (minusgood)
        gluc = (glucminus+gluc)/2;
    else
        error('no good solution')
    end       
    
    ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
    ksx_int = dx * (sum(ksx)-(ksx(1)+ksx(end))/2);
    if (opt.highks)
         M  = avgmitodenh  * (ksx/kwh) ./ ((ksx_int/kwh)) + eps;
    else
        M  = avgmitodenh  * (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    end
    normdg = norm(glucold - gluc);
    step = step+1;
    
    if (step>opt.nstep)
        disp('Failed to converge')
        ftc = 1;
        return
    end
    % plotting options
    if (opt.dodisplay && mod(step,opt.showevery)==0)
        
        plot(xpos,gluc_init/(opt.L)^3,'k--')
        hold all
        plot(xpos,gluc/(opt.L)^3,'b.-')
        plot(xpos,C_out,'g--')
        plot(xpos,M/(opt.L)^2 ,'r.-')
        title(sprintf('Step %d, normdg/gluc: %f', step, normdg/mean(gluc)))
        hold off
        drawnow
        
    end
end


        if (opt.dodisplay)
        plot(xpos,gluc_init/(opt.L)^3,'k--')
        hold all
        plot(xpos,gluc/(opt.L)^3,'b.-')
        plot(xpos,C_out/(opt.L)^3,'g--')
        plot(xpos,M/(opt.L)^3,'r.-')
        title(sprintf('Step %d, normdg/gluc: %f', step, normdg/mean(gluc)))
        hold off
        drawnow
        end
        
    
end
    


