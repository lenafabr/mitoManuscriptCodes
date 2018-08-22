%code to produce discrete simulation plots 
function  [gluc, mitopos, mitostate, opt] = rundiscretesims(options)
% NOTE: this now assumes concentration units are provided in #/length^3

%% set up default simulation parameters
opt = struct();

opt.kg = 1; % rate of glucose consumption
opt.c0 = 1; % fixed glucose concentration
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.r = 0.4; % domain radius 
opt.vel = 1; % mitochondria velocity
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 200; % rate of stopping is ks*[gluc] 
opt.Km = 0.1; % Michaelis constant for glucose consumption

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

% fix permanent glucose distribution; do not evolve it
opt.fixgluc = [];

opt.nmito = 75; % number of mitochondria, updated from Gulcin's paper
opt.gpts = 500; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-3; % time-step
opt.nstep = 1e7; % number of steps to run

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% prob mitochondria start in walking state
% default is use equilibrium probability at c0 conc
opt.pstartwalk = opt.kw/(opt.kw + opt.ks*opt.c0);

% fix glucose concentrations
opt.glucfix = 0;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 5e3;
opt.showmito = 1;

opt.restart = 1; % flag to enable continuing previous sims

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

% boundary conditions on the far size
% positive = fixed concentration at the boundary
% negative = reflecting boundary
if(~isfield(opt,'cend'))
    opt.cend = opt.c0;
end

% set up dimensionless parameters
%dimensionless parameters are L(length unit),L^2/D(time units)
tscale = (opt.L)^2 / opt.D;
lscale = opt.L;

Lh = opt.L/lscale;
rh = opt.r/lscale;

kwh = opt.kw*tscale;
ksh = opt.ks*tscale;
Dh = opt.D*tscale/lscale^2;
kgh = opt.kg*tscale;
c0h = opt.c0*lscale^3;
msizeh = opt.msize/lscale;
velh = opt.vel*tscale/lscale;
cendh = opt.cend*lscale^3;
startpos = opt.startpos/lscale;
Kmh = opt.Km*lscale^3;

delth = opt.delt/tscale;
% spatial resolution
dx = Lh/(opt.gpts - 1);

%%
% initialize variables


% spatial positions at which glucose is evaluated
xpos = linspace(0,Lh,opt.gpts)';

if (opt.restart)
    % mitochondria center positions
    if (length(opt.startpos) == opt.nmito)
        % start with predefined positions
        mitopos = startpos;
    elseif (length(opt.startpos)==2)
        % start evenly split between specific positions
        u = rand(opt.nmito,1);
        mitopos = zeros(opt.nmito,1);
        mitopos(u<0.5) = startpos(1);
        mitopos(u>=0.5) = startpos(2);
    elseif opt.startpos>0
        % start at specific position
        mitopos = startpos*ones(opt.nmito,1);
    else
        % start uniformly distributed btwn msizeh and Lh-msizeh
        mitopos = rand(opt.nmito,1)*(Lh-2*msizeh)+msizeh;
    end
    
    % glucose concentration
    % start with linear conc profile from c0 to cend
    if opt.cend<0
        error('reflecting boundary not yet implemented')
    else
        if (isempty(opt.startgluc))
            gluc = linspace(c0h,cendh,opt.gpts)';
        else
            if (length(opt.startgluc) ~= opt.gpts); error('starting distrib has wrong size'); end
            gluc = opt.startgluc;
        end
    end
    
    
    % state of the mitochondria (walking or not)
    % each mitochondria starts walking with probability p
    u = rand(opt.nmito,1);
    mitostate = u<=opt.pstartwalk;
    
    % set initial direction of walking (randomly)
    u = rand(opt.nmito,1);
    mitostate = ((u<=0.5)*2-1).*mitostate;
end

if (~isempty(opt.fixgluc))
    if (length(opt.fixgluc) ~= length(xpos)); error('wrong length of fixgluc'); end
    gluc = opt.fixgluc;
end

mitopos0 = mitopos;
%% evolve the system over time
d2g = zeros(opt.gpts,1);
dtg = d2g;

% probability of starting on each timestep
pstart = 1 - exp(-kwh*delth);

if (opt.restart); curtime = 0; end

for step = 1:opt.nstep
    
    if (any(mitopos>Lh-0.5*msizeh) || any(mitopos<0.5*msizeh))
        error('bad mito positions')
    end
    
    if (isempty(opt.fixgluc))
        % ---------
        % evolve forward the glucose concentration by 1 time step
        % ----------
        if (~opt.glucfix)
            % second derivative of gluc
            d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2;
            % time derivative due to diffusion
            dtg = Dh*d2g;
        end
        % time derivative due to glucose consumption
        % estimate consumption rate at all the mito positions
        %glucmito = interp1(xpos,gluc,mitopos);
        %consum = kgh*Kmh*glucmito./(Kmh+glucmito);
        
        glucmito = zeros(opt.nmito,1);
        for mc = 1:opt.nmito            

            
            ind1 = ceil((mitopos(mc)-1/2*msizeh)/dx)+1;
            ind2 = floor((mitopos(mc)+1/2*msizeh)/dx)+1;
            glucmito(mc) = (gluc(ind1)+gluc(ind2))/2;
            % dimensionless consumption rate of 1
            % note: overlapping mitochondria will consume twice as fast
            if (ind2>=ind1)
                dtg(ind1:ind2) = dtg(ind1:ind2)-kgh/(pi*rh^2*msizeh)*gluc(ind1:ind2)./(Kmh+gluc(ind1:ind2));
            end
%             %----------
        end
        
        if (~opt.glucfix)
            % fix boundary conditions
            dtg(1) = 0;
            if (opt.cend>0)
                % fixed conc at far end
                dtg(end) = 0;
            else
                % reflecting boundary at far end
                error('reflecting boundary not yet implemented')
            end
            
            % evolve forward
            gluc = gluc+dtg*delth;
        end
        
        if (any(gluc<-1e-3))
            error('negative concentrations!')
        end
    end
    
    % move the walking mitochondria
    walkind = find(mitostate);
    stopind = find(~mitostate);
    mitopos(walkind) = mitopos(walkind) + velh*delth*mitostate(walkind);
    
    %% reflect mitochondria back if hitting the boundary
    for mc = walkind'
        if (mitopos(mc)<0.5*msizeh)
            mitostate(mc)=1;
            mitopos(mc) = 0.5*msizeh+(0.5*msizeh-mitopos(mc));
        elseif mitopos(mc)>Lh-0.5*msizeh
            mitostate(mc) = -1;
            mitopos(mc) = Lh-0.5*msizeh - (mitopos(mc)-Lh+0.5*msizeh);
        end
    end
    
    %%
    % decide which mitochondria stop
    % glucose concentrations at mitochondria positions
    % use glucose at discrete points around mito position
    %glucmito = interp1(xpos,gluc,mitopos(walkind));
    %stoprate = ksh*Kmh*glucmito./(Kmh+glucmito);
    
    stoprate = ksh*glucmito(walkind)./(Kmh+glucmito(walkind));
    pstop = 1-exp(-stoprate*delth);
    u = rand(length(walkind),1);
    mitostate(walkind) = mitostate(walkind).*(1 - (u<=pstop));
    
    % decide which mitochondria start after this step
    u = rand(length(stopind),1);
    restartind = find(u<=pstart);
    mitostate(stopind(restartind)) = (rand(length(restartind),1)<=0.5)*2 - 1; % start in random direction
    
    if (opt.dodisplay && mod(step,opt.showevery)==0)
        % plot glucose concentration
        plot(xpos,gluc,'.-')
        ylim([0,c0h])
        hold on;
        % plot mitochondria positions
        if (opt.showmito)
            hold all
            ymin = min(gluc); ymax = max(gluc);
            for mc = 1:opt.nmito
                if (mitostate(mc)==0)
                    
                    plot([mitopos(mc)], [0.5*c0h],'or','LineWidth',2)
                else
                    
                    plot([mitopos(mc)], [0.5*c0h],'ob','LineWidth',2)
                end
                set(gca,'FontSize',12)
                
            end
        end
        mval = (6*var(mitopos)/Lh^2) - 0.5;
        title(sprintf('Step %d: variance metric = %0.3f',step,mval))
        hold off
 
        drawnow
    end
    
    curtime = curtime +delth;
end
end
