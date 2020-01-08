function [P1,sett,model] = GetModel(model_num,P,ix,dir_res,opt)
% Gey models

numsubj = opt.numsubj;
run3d   = opt.run3d;
ax2d    = opt.ax2d;

model = [];
P1    = {};
sett  = struct;

if model_num == 0, fprintf('=============\nMODEL 0\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 0 | Testing
%------------------

% Define training population
ix_pop = ix.IXI;
N      = 8; % Set maximum number of subjects
P1     = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N;
    P1{p}{4} = 1;
end

% ix_pop = [ix.BALGRIST ix.IXI ix.MICCAI2012 ix.MRBRAINS18]; % Set training populations to use
% N      = 2; % Set maximum number of subjects
% 
% % Define training population
% P1 = P(ix_pop);
% for p=1:numel(P1)
%     P1{p}{2} = {'T1'};
%     P1{p}{3} = N;
%     P1{p}{4} = 1;
% end

% ix_pop = ix.IXIRC;
% N      = 4; % Set maximum number of subjects
% 
% % Define training population
% P1       = P(ix_pop);
% P1{1}{3} = N;

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','normalised','intensity'};
sett.write.dir_res      = fullfile(dir_res,'results/model-0');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.model.mg_ix        = [1 1 1 1 2 3 4 5 6];
sett.labels.use         = true; 
sett.model.K            = 6;  
sett.show.mx_subjects   = 8;
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 1 || model_num == 2
%%%%%%%%%%%%%%%%%%%
% Model 1 and 2
%------------------

% Set training populations to use
ixs    = [ix.IXI ix.MRBRAINS18 ix.BALGRIST ix.DELIRIUM ix.MICCAI2012 ix.ATLAS];
N      = [80 4 12 80 20 80]; % Set maximum number of subjects
N      = min(N,numsubj);
int_ix = [1 3 2 4 3 3];

% Label information
igm = 8; iwm = 6; icsf = 9; iven = 9; k1 = 10; % OBS: change csf 10 -> 9 for 3D
cm_map = {{}, {igm,iwm,icsf,iven,[igm iwm],setdiff(1:k1,[iven])}, {iwm,[]}, {}, {igm,iwm,iven,setdiff(1:k1,[iven])}, {}};

% Define training population
P1 = P(ixs);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_ix(p);
    P1{p}{5} = cm_map{p};
end
P1{2}{2} = {'T1'};

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','intensity'};
sett.model.init_mu_dm   = 32;
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.mu           = [true true];
sett.nit.init           = 6;
sett.show.mx_subjects   = 4;
end

if model_num == 2, fprintf('=============\nMODEL 1\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 1 | Labels are used (K1=10), trying to get nice GM, WM and CSF
%------------------

% Settings
sett.write.dir_res = fullfile(dir_res,'results/model-1');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use    = true; 
sett.model.K       = 9;  
sett.model.mg_ix   = [1 1 2 2 3 4 5 6 7 8 9 10];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 3, fprintf('=============\nMODEL 2\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 2 | Labels are not used (K1=12), unsupervised for better normalisation
%------------------

% Settings
sett.write.dir_res = fullfile(dir_res,'results/model-2');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use    = false; 
sett.model.K       = 11;  
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 4, fprintf('=============\nMODEL 3\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 3 | Labels are used (K1=7, mg_ix=2), trying to get nice GM, WM and CSF
%------------------

% Set training populations to use
ixs    = [ix.IXI ix.MRBRAINS18 ix.BALGRIST ix.DELIRIUM ix.MICCAI2012 ix.ATLAS];
N      = [80 4 12 80 20 80]; % Set maximum number of subjects
N      = min(N,numsubj);
int_ix = [1 5 2 4 3 3];

% Label information
igm = 6; iwm = 5; icsf = 3; iven = 3; k1 = 8;
cm_map = {{}, {igm,iwm,icsf,iven,[igm iwm],setdiff(1:k1,[iven])}, {iwm,[]}, {}, {igm,iwm,iven,setdiff(1:k1,[iven])}, {}};

% Define training population
P1 = P(ixs);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_ix(p);
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-3');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = 7;  
sett.model.mg_ix        = 2;
sett.model.ix_init_pop  = 5;
sett.write.mu           = [true true];
sett.nit.init           = 6;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
end

if model_num == 5, fprintf('=============\nMODEL 4\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 4 | Fit only T1 (use to InitGMM, K1=12), unsupervised
%------------------

% Set training populations to use
ixs    = [ix.MRBRAINS18  ix.IXI ix.BALGRIST ix.MICCAI2012];
N      = [4 80 12 80 20]; % Set maximum number of subjects
N      = min(N,numsubj);
int_ix = [1 1 1 1];
cm_map = {{},{},{},{}};

% Define training population
P1 = P(ixs);
for p=1:numel(P1)    
    P1{p}{2} = {'T1'};
    P1{p}{3} = N(p);
    P1{p}{4} = int_ix(p);
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-4');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = false; 
sett.model.K            = 11;  
sett.model.mg_ix        = 1;
sett.model.ix_init_pop  = 1;
sett.write.mu           = [true true];
sett.nit.init           = 6;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 6, fprintf('=============\nMODEL 5\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 5 | Fit a learned model to new subjects
%------------------

% Define training population
ix_pop   = ix.MICCAI2012; % Set training populations to use
P1       = P(ix_pop);
P1{1}{2} = {'T1'};
P1{1}{3} = {21:30};
P1{1}{4} = 3;

% Settings
sett                   = struct;
sett.show.figs         = {'model','segmentations','normalised','parameters','intensity'};
sett.write.dir_res     = fullfile(dir_res,'results/model-5');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
sett.write.tc          = false(1,3);
sett.write.im          = false(1,4);
sett.write.df          = [true true];
sett.write.labels      = [true false];
sett.do.infer          = true;
sett.model.init_mu_dm  = 16;  
sett.var.v_settings    = [1e-4 0 0.2 0.05 0.2]*2^3;
sett.show.mx_subjects  = 4;
sett.clean_z.mrf       = 2;
sett.clean_z.gwc_tix   = struct('gm',[8 10],'wm',[7],'csf',[11]);

% Path to a model
% model_num = 1;
% pth_model = fullfile(dir_res,['results/model-' num2str(model_num)]);
% if ~run_3D, pth_model = [pth_model '-2D-' ax_2D]; end
% pth_model = fullfile(pth_model,'model_spm_mb.mat');
pth_model = '/scratch/Results/diffeo-segment/k11/model_spm_mb.mat';

% Load model
load(pth_model) % loads variable 'model'
% model = rmfield(model,'shape');  % uncommented -> no mu
% model = rmfield(model,'appear'); % uncommented -> no pr
end
end
%==========================================================================