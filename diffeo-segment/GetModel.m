function [P1,sett,model] = GetModel(model_num,P,ix,dir_res,opt)
% Get models

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
ix_pop = [ix.MICCAI2012];
N      = [8 100 100 4 12 80];
N      = min(N,numsubj);
int_ix = [1 2];

% Number of template classes
K = 5;

% Label information
% icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5; k1 = K + 1;
% cm_map = {{icgm,isgm,[isgm iwm],iwm,iven, setdiff(1:k1,[icgm isgm iven])}, ...
%           {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:k1,[icgm isgm iven])}, ...
%           {iwm,[]}};
cm_map = cell(1,numel(ix_pop));

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = int_ix(p);      
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','parameters','intensity'};
sett.write.dir_res      = fullfile(dir_res,'results/model-0');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.model.mg_ix        = 1;
sett.labels.use         = false; 
sett.model.K            = K;  
sett.show.mx_subjects   = 2;
sett.model.init_mu_dm   = 16;
sett.nit.init           = 6;
sett.write.mu           = [true true];
sett.var.v_settings     = [1e-4 0 0.2 0.05 0.2]*2;     
sett.optim.scal_q       = 1;
sett.optim.scal_v       = 0.5;
sett.optim.scal_bf      = 1;
sett.gen.samp           = 4;
sett.gen.samp_mx        = 1;
sett.do.updt_aff        = true;
sett.do.updt_bf         = true;
sett.do.updt_int        = true;
sett.do.updt_vel        = true;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

if model_num == 1, fprintf('=============\nMODEL 1\n=============\n\n');
%%%%%%%%%%%%%%%%%%%
% Model 1 | Template prop MRBRAINS18, K1=12, get CGM, SGM, WM, CSF, VEN
%------------------

% Define training population
ix_pop  = [ix.MICCAI2012 ix.MRBRAINS18 ix.BALGRIST  ix.IXI ix.ROB];
N       = [100 100 100 100 100]; % Set maximum number of subjects
N       = min(N,numsubj); 
int_pop = [1 1 2 3 4];

% Number of template classes
K = 11;

% Label information
icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5; k1 = K + 1;
cm_map = {{icgm,isgm,[isgm iwm],iwm,icsf,iven, setdiff(1:k1,[icgm])}, ...
          {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:k1,[icgm])}, ...
          {[isgm iwm],[]},{},{}};
% cm_map = cell(1,numel(ix_pop));

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_pop(p);  
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.model.init_mu_dm   = 32;  
sett.write.intermediate = false;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-1');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = K;  
sett.model.mg_ix        = 1;
sett.model.ix_init_pop  = 1;   
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
N      = min(N,repmat(numsubj,[1 numel(ixs)]));

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
% Model 3 | Template prop BALGRIST, K1=8, get GM, WM, CSF
%------------------

% Define training population
ix_pop  = [ix.MRBRAINS18 ix.MICCAI2012 ix.BALGRIST ix.IXI ix.ROB];
N       = 100*ones(1,numel(ix_pop)); % Set maximum number of subjects
N       = min(N,numsubj); 
int_pop = [1 2 3 4 5];

% Number of template classes
K  = 7;
K1 = K + 1;

% Label information
igm = 1; iwm = 2; icsf = 3;
cm_map = {{igm,igm,[igm iwm],iwm,[igm iwm],icsf,icsf,1:K1}, ...
          {igm,igm,[igm iwm],iwm,icsf,icsf, 1:K1}, ... 
          {[igm iwm],[]},...                   
          {},{}};
% cm_map = cell(1,numel(ix_pop));

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_pop(p);  
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = false;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,'results/model-3');
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = K;  
sett.model.ix_init_pop  = 1;   
sett.show.mx_subjects   = 2;
sett.write.mu           = [true true];
sett.nit.init           = 8;
sett.var.v_settings     = [1e-4 0 0.2 0.05 0.2]*2^2;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
end
end
%==========================================================================