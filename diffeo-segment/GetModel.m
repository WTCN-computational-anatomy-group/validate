function [P1,sett,model] = GetModel(model_num,P,ix,dir_res,opt)
% MODELS:
% 0. Testing
% 1. Template prop MICCAI2012+MRBRAINS18, (K1=12), get CGM, SGM, WM, CSF, VEN
% 2. Fit T1w (K1=12)
% 3. MRI+CT model (K1=12)
% 4. CROMIS (K1=12)
% 5. MICCAI2020 (IBSR18)
% 6. MICCAI2020 (LPBA40)
% 7. MICCAI2012
% 8. CROMISPETTERI
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

numsubj = opt.numsubj;
run3d   = opt.run3d;
ax2d    = opt.ax2d;
nw      = opt.nw;
showfig = opt.show;

model = [];
P1    = {};
sett  = struct;

%%%%%%%%%%%%%%%%%%%
% Model 0 | Testing
%------------------
if model_num == 0, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop = [ix.IXI];
N      = numsubj*ones(1,numel(ix_pop));

% Number of template classes
K = 9; K1 = K + 1;

% Label information
% icgm   = 1; isgm = 2; iwm = 3; icsf = 4; iven = 5;
% cm_map = {{icgm,isgm,[isgm iwm],iwm,iven, setdiff(1:K1,[icgm isgm iven])}, ...
%           {icgm,isgm,[isgm iwm],iwm,[isgm iwm],icsf,iven,setdiff(1:K1,[icgm isgm iven])}, ...
%           {iwm,[]}};

cm_map = cell(1,numel(ix_pop));
int_ix = 1:numel(ix_pop);
P1     = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = int_ix(p);      
    P1{p}{5} = cm_map{p};
end
P1{1}{2} = P1{1}{2}(1:3);

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations'}; end
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
sett.labels.use         = false; 
sett.model.K            = K;  
sett.write.scal_mom     = true;
sett.do.updt_bf         = true;
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
K = 11; K1 = K + 1;

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
% P1{4}{2} = {'T1','T2','PD'};

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations','init','intensity'}; end
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory 
sett.write.mu           = [true true];
sett.write.workspace    = true;
sett.write.df           = true;
sett.labels.use         = true; 
sett.model.K            = K;  
sett.model.mg_ix        = 1;
sett.model.ix_init_pop  = 1;   
sett.show.mx_subjects   = 8;
sett.model.crop_mu      = false;
end
%%%%%%%%%%%%%%%%%%%
% Model 2 | Fit T1w (K1=12)
%------------------
if model_num == 2, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Set training populations to use
ixs = [ix.IXI ix.BALGRIST ix.MADRID ix.MICCAI2012 ix.MRBRAINS18];
N   = [200 Inf Inf Inf Inf];
N   = min(N,numsubj);

% Number of template classes
K = 11; K1 = K + 1;

% Define training population
P1 = P(ixs);
for p=1:numel(P1)    
    P1{p}{2} = {'T1'};
    P1{p}{3} = N(p);
    P1{p}{4} = 1;
end

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations','init','intensity'}; end
sett.show.mx_subjects   = 8;
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.workspace    = true;
sett.write.mu           = [true true];
sett.write.df           = true;
sett.labels.use         = false; 
sett.model.K            = K;  
sett.model.crop_mu      = false;
sett.model.init_mu_dm   = 8;
sett.nit.init           = 32;
sett.var.v_settings     = [0 0 0.2 0.05 0.2]*2;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end
%%%%%%%%%%%%%%%%%%%
% Model 3 | MRI+CT model (K1=12)
%------------------
if model_num == 3, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop  = [ix.BALGRIST ix.IXI ix.MICCAI2012 ix.MRBRAINS18 ix.CTHEALTHY ix.ROB ix.MPMCOMPLIANT]; 
N       = [Inf 300 Inf Inf Inf Inf Inf]; 
N       = min(N,numsubj); 
int_pop = [1 2 3 3 4 4 5]; 

% Number of template classes
K = 12; K1 = K + 1;

igm = 13; iwm = 12; icsf = 8; ibg = 2; ibn = 5;
cm_map = {{[igm iwm],1:K1},... 1.spn 
          {}, ...
          {igm,igm,[igm iwm],iwm,icsf,icsf,1:K1}, ... % 1.cgm,2.sgm,3.spn,4.wm,5,csf,6.ven  
          {igm,igm,[igm iwm],iwm,[igm iwm],icsf,icsf,1:K1}, ...  % 1.cgm,2.sgm,3.spn,4.wm,5.cer,6.csf,7.ven                                     
          {ibg,ibn,setdiff(1:K1,ibn)}, ...  % 1.bg,2.bone
          {},{}};
          
P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p);
    P1{p}{4} = int_pop(p);  
    P1{p}{5} = cm_map{p};
end
if numel(P1) > 1, P1{2}{2} = {'T1','T2','PD'}; end

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations','init','intensity','parameters'}; end
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df           = false;
sett.write.mu           = [true true];
sett.labels.use         = true; 
sett.labels.use_initgmm = false;
sett.model.K            = K;  
sett.model.ix_init_pop  = 1;
sett.model.crop_mu      = false;
sett.nit.init           = 32;
sett.nit.zm             = 4;
sett.model.mg_ix        = 1;
sett.model.init_mu_dm   = 32;
sett.var.v_settings     = [0 0 0.2 0.05 0.2]*2;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end
%%%%%%%%%%%%%%%%%%%
% Model 4 | CROMIS (K1=12)
%------------------
if model_num == 4, fprintf('=============\nMODEL %i\n=============\n\n',model_num);

% Define training population
ix_pop = [ix.CROMIS];
N      = Inf*ones(1,numel(ix_pop)); % Set maximum number of subjects
N      = min(N,numsubj);

% Number of template classes
K = 11; K1 = K + 1;
 
cm_map = cell(1,numel(ix_pop));
% ibg = 1; ibn = K1;
% cm_map = {{ibg,ibn,setdiff(1:K1,ibn)}};

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = 1;      
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations','init','intensity'}; end
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df           = true;
sett.write.mu           = [true true];
sett.write.tc           = [true false false];  % native, warped, warped-mod
sett.model.K            = K;  
sett.show.mx_subjects   = 8;
sett.nit.init           = 32;
sett.nit.zm             = 3;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end
%%%%%%%%%%%%%%%%%%%
% Model 5 | MICCAI2020 (IBSR18)
%------------------
if model_num == 5, fprintf('=============\nMODEL %i\n=============\n\n',model_num);
% Set training populations to use
ixs = ix.IBSR18;
N   = 18;
N   = min(N,numsubj);

% Number of template classes
K = 11; K1 = K + 1;

% Define training population
P1       = P(ixs);   
P1{1}{3} = N;

% Settings
sett                  = struct;
if showfig, sett.show.figs = {'model','segmentations'}; end
sett.show.mx_subjects = 8;
sett.gen.num_workers  = nw;
sett.write.dir_res    = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df         = true;
sett.write.mu         = [true true];
sett.model.K          = K;
sett.model.init_mu_dm = 8;
sett.var.v_settings   = [0 0 0.2 0.05 0.2]*2;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory    
end
%%%%%%%%%%%%%%%%%%%
% Model 6 | MICCAI2020 (LPBA40)
%------------------
if model_num == 6, fprintf('=============\nMODEL %i\n=============\n\n',model_num);
% Set training populations to use
ixs = ix.LPBA40;
N   = 40;
N   = min(N,numsubj);

% Number of template classes
K = 11; K1 = K + 1;

% Define training population
P1       = P(ixs);   
P1{1}{3} = N;

% Settings
sett                  = struct;
if showfig, sett.show.figs = {'model','segmentations'}; end
sett.show.mx_subjects = 8;
sett.gen.num_workers  = nw;
sett.write.dir_res    = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df         = true;
sett.write.mu         = [true true];
sett.model.K          = K;
sett.model.init_mu_dm = 8;
sett.var.v_settings   = [0 0 0.2 0.05 0.2]*2;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory    
end
%%%%%%%%%%%%%%%%%%%
% Model 7 | MICCAI2012
%------------------
if model_num == 7, fprintf('=============\nMODEL %i\n=============\n\n',model_num);
% Set training populations to use
ixs = ix.MICCAI2012;
N   = 35;
N   = min(N,numsubj);

% Number of template classes
K = 11; K1 = K + 1;

% Define training population
P1       = P(ixs);   
P1{1}{3} = N;

% Settings
sett                  = struct;
if showfig, sett.show.figs = {'model','segmentations'}; end
sett.show.mx_subjects = 8;
sett.gen.num_workers  = nw;
sett.write.dir_res    = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df         = true;
sett.write.mu         = [true true];
sett.model.K          = K;
sett.model.init_mu_dm = 8;
sett.var.v_settings   = [0 0 0.2 0.05 0.2]*1.5;
sett.nit.init         = 128;
sett.nit.zm           = 4;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory    
end
%%%%%%%%%%%%%%%%%%%
% Model 8 | CROMISPETTERI
%------------------
if model_num == 8, fprintf('=============\nMODEL %i\n=============\n\n',model_num);
% Define training population
ix_pop = [ix.CROMISPETTERI];
N      = Inf*ones(1,numel(ix_pop)); % Set maximum number of subjects
N      = min(N,numsubj);

% Number of template classes
K = 10; K1 = K + 1;
 
ibrain = 1:5;
iven = 6;
iles = 7;
ical = 8;
ibn = 9;
irest = 10;
ibg = 11;     
cm_map = {{ibg, ibn, iven, iles, ical, ibrain, irest, ...
           setdiff(1:K1,[iven, ical, iles])}};

P1 = P(ix_pop);
for p=1:numel(P1)
    P1{p}{3} = N(p); 
    P1{p}{4} = 1;      
    P1{p}{5} = cm_map{p};
end

% Settings
sett                    = struct;
if showfig, sett.show.figs = {'model','segmentations','intensity'}; end
sett.gen.num_workers    = nw;
sett.write.dir_res      = fullfile(dir_res,['results/model-' num2str(model_num)]);
if ~run3d, sett.write.dir_res = [sett.write.dir_res '-2D-' ax2d]; end
sett.write.df           = false;
sett.write.mu           = [true true];
sett.write.tc           = [false false false];  % native, warped, warped-mod
sett.model.K            = K;  
sett.nit.zm             = 3;
sett.model.mg_ix        = [1 2 3 4 5 6 7 8 9 10 10 11 11];
sett.model.vx           = 1;
sett.labels.use         = true; 
sett.labels.use_initgmm = true;
sett.var.v_settings     = [0 0 0.2 0.05 0.2]*4;
sett.nit.init           = 16;
sett.model.init_mu_dm   = 16;
sett.do.updt_bf         = true;
sett.labels.w           = 0.9999;
sett.model.crop_mu      = false;
if exist(sett.write.dir_res,'dir') == 7, rmdir(sett.write.dir_res,'s'); end % clear results directory
end
%==========================================================================
