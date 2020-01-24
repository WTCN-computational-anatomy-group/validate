function [P1,sett,model] = GetModel(model_num,P,ix,dir_res,opt)
% MODELS:
% 0. Testing
% 1. Template prop MICCAI2012+MRBRAINS18, (K1=12), get CGM, SGM, WM, CSF, VEN
% 2. Fit only T1 (K1=12)
% 3. Template prop MRBRAINS18, (K1=8), get GM, WM, CSF
% 4. Fit IXI (K1=14)
% 5. CROMIS (K1=12)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

numsubj = opt.numsubj;
run3d   = opt.run3d;
ax2d    = opt.ax2d;

model = [];
P1    = {};
sett  = struct;

%%%%%%%%%%%%%%%%%%%
% Model 0 | Testing
%------------------
if model_num == 0, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop = [ix.MICCAI2012];
N      = 100*ones(1,numel(ix_pop));
N      = min(N,numsubj);
int_ix = [1 2];

% Number of template classes
K  = 5;
K1 = K + 1;

% Label information
% icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5;
% cm_map = {{icgm,isgm,[isgm iwm],iwm,iven, setdiff(1:K1,[icgm isgm iven])}, ...
%           {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:K1,[icgm isgm iven])}, ...
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
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.model.mg_ix        = 1;
sett.labels.use         = false; 
sett.model.K            = K;  
sett.show.mx_subjects   = 2;
sett.write.mu           = [true true];
sett.var.v_settings     = [1e-4 0 0.2 0.05 0.2]*4;     
sett.do.updt_vel        = true;
sett.do.updt_aff        = true;
sett.do.updt_int        = true;
sett.do.updt_template   = true;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

%%%%%%%%%%%%%%%%%%%
% Model 1 | Template prop MICCAI2012+MRBRAINS18, (K1=12), get CGM, SGM, WM, CSF, VEN
%------------------
if model_num == 1, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop  = [ix.MICCAI2012 ix.MRBRAINS18 ix.BALGRIST  ix.IXI ix.ROB];
N       = 100*ones(1,numel(ix_pop)); % Set maximum number of subjects
N       = min(N,numsubj); 
int_pop = [1 1 2 3 4];

% Number of template classes
K  = 11;
K1 = K + 1;

% Label information
icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5;
cm_map = {{icgm,isgm,[isgm iwm],iwm,icsf,iven, setdiff(1:K1,[icgm])}, ...
          {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:K1,[icgm])}, ...
          {[isgm iwm],[]},{},{}};
% cm_map = cell(1,numel(ix_pop));

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_pop(p);  
    P1{p}{5} = cm_map{p};
end
P1{4}{2} = {'T1','T2','PD'};

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.model.init_mu_dm   = 32;  
sett.write.intermediate = false;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = K;  
sett.model.mg_ix        = 1;
sett.model.ix_init_pop  = 1;   
sett.show.mx_subjects   = 2;
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
end

%%%%%%%%%%%%%%%%%%%
% Model 2 | Fit only T1 (K1=12)
%------------------
if model_num == 2, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Set training populations to use
ixs    = [ix.IXI  ix.MICCAI2012 ix.BALGRIST];
N      = [150 20 10]; % Set maximum number of subjects
N      = min(N,numsubj);

% Number of template classes
K  = 11;
K1 = K + 1;

% Define training population
P1 = P(ixs);
for p=1:numel(P1)    
    P1{p}{2} = {'T1'};
    P1{p}{3} = N(p);
    P1{p}{4} = p;
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = false; 
sett.model.K            = K;  
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

%%%%%%%%%%%%%%%%%%%
% Model 3 | Template prop MRBRAINS18, (K1=8), get GM, WM, CSF
%------------------
if model_num == 3, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

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
P1{4}{2} = {'T1','T2','PD'};

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = false;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = true; 
sett.model.K            = K;  
sett.model.ix_init_pop  = 1;   
sett.show.mx_subjects   = 2;
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
end

%%%%%%%%%%%%%%%%%%%
% Model 4 | Fit IXI (K1=14)
%------------------
if model_num == 4, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

N = 100; % Set maximum number of subjects
N = min(N,numsubj);

% Number of template classes
K  = 13;
K1 = K + 1;

% Define training population
P1       = P(ix.IXI);  
% P1{1}{2} = {'T1','T2','PD'}; % uncomment to not include MRA
P1{1}{3} = N;

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.model.init_mu_dm   = 16;  
sett.write.intermediate = true;
sett.write.clean_vel    = false;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.labels.use         = false; 
sett.model.K            = K;  
sett.write.mu           = [true true];
sett.model.crop_mu      = false;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

%%%%%%%%%%%%%%%%%%%
% Model 5 | CROMIS (K1=12)
%------------------
if model_num == 5, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop = [ix.CROMISLABELS ix.CROMIS];
N      = Inf*ones(1,numel(ix_pop)); % Set maximum number of subjects
N      = min(N,numsubj);

% Number of template classes
K  = 11;
K1 = K + 1;
 
% Label information
iles   = 1;
cm_map = {{iles, setdiff(1:K1,iles)}, {}};

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = 1;      
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
sett.show.figs          = {'model','segmentations','InitGMM','intensity'};
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.intermediate = true;
sett.model.mg_ix        = 1;
sett.labels.use         = true; 
sett.model.K            = K;  
sett.show.mx_subjects   = 4;
sett.write.mu           = [true true];
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end

end
%==========================================================================