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
ix_pop = [ix.MICCAI2012 ix.MRBRAINS18];
N      = [100 100 100 4 12 80];
N      = min(N,numsubj);
int_ix = [1 1 2];

% Label information
icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5; k1 = 10;
cm_map = {{icgm,isgm,[isgm iwm],iwm,iven, setdiff(1:k1,[icgm isgm iven])}, ...
          {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:k1,[icgm isgm iven])}, ...
          {iwm,[]}};

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = int_ix(p);      
    P1{p}{5} = cm_map{p};
end
% P1{3}{2} = {'T1'};

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.write.dir_res      = fullfile(dir_res,'results/model-0');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.model.mg_ix        = 1;
sett.labels.use         = false; 
sett.model.K            = 13;  
sett.show.mx_subjects   = 6;
sett.model.init_mu_dm   = 16;
sett.nit.init           = 6;
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 1, fprintf('=============\nMODEL 1\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 1 | Fit K1=8 supervised, with MRI and CT, try to get nice tissues
%------------------

% Define training population
ix_pop  = [ix.MICCAI2012 ix.IXI ix.BALGRIST ix.ROB];
N       = [100 100 100 100]; % Set maximum number of subjects
N       = min(N,numsubj); 
int_pop = [1 2 3 4];

% Label information
icgm   = 4; isgm = 6; iwm = 3; icsf = 5; k1 = 8;
cm_map = {{icgm,isgm,[isgm iwm],iwm,icsf, setdiff(1:k1,[icgm isgm])}, {}, {iwm,[]}, {}};

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_pop(p);  
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-1');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = 7;  
sett.model.mg_ix        = [1 1 1 2 3 4 5 5 6 7 8];
sett.show.mx_subjects   = 2;
sett.write.mu           = [true true];
sett.nit.init           = 6;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
end

if model_num == 2, fprintf('=============\nMODEL 2\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 2 | Fit only T1 (use to InitGMM, K1=14), unsupervised
%------------------

% Set training populations to use
ixs    = [ix.IXI  ix.MICCAI2012 ix.BALGRIST];
N      = [150 20 10]; % Set maximum number of subjects
N      = min(N,numsubj);

% Define training population
P1 = P(ixs);
for p=1:numel(P1)    
    P1{p}{2} = {'T1'};
    P1{p}{3} = N(p);
    P1{p}{4} = 1;
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-2');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = false; 
sett.model.K            = 13;  
sett.write.mu           = [true true];
sett.nit.init           = 6;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 3, fprintf('=============\nMODEL 3\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 3 | Fit a learned model to new subjects
%------------------

% Define training population
ix_pop   = ix.MICCAI2012; % Set training populations to use
P1       = P(ix_pop);
P1{1}{2} = {'T1'};
P1{1}{3} = {21};
P1{1}{4} = 3;

% Settings
sett                   = struct;
sett.show.figs         = {'model','segmentations','normalised','parameters','intensity'};
sett.write.dir_res     = fullfile(dir_res,'results/model-3');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
sett.write.tc          = [true true false];
sett.write.im          = [false false false true];
sett.write.df          = [true true];
sett.write.labels      = [true false];
sett.do.infer          = true;
sett.model.init_mu_dm  = 16;  
sett.var.v_settings    = [1e-4 0 0.2 0.05 0.2]*2^2;
sett.show.mx_subjects  = 4;
% sett.clean_z.mrf       = 2;
% sett.clean_z.gwc_tix   = struct('gm',[8 10],'wm',[7],'csf',[11]);
sett.write.vel         = true;
sett.write.affine      = true;

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