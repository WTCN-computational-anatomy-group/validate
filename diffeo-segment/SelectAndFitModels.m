function SelectAndFitModels(opt)
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
%
% SelectAndFitModels(struct('user','mbrud-fil','models',0,'run3d',true,'numsubj',Inf,'nw',0))
%
% -------------------------------------------------------------------------
% TODO
% -------------------------------------------------------------------------
% [ ] profile Register
% [ ] tissue prop
% [ ] Use binning uncertainity instead of jitter? Worth it? Slower right?
%
% -------------------------------------------------------------------------
% QUESTIONS
% -------------------------------------------------------------------------
% [ ]
%
% -------------------------------------------------------------------------
% EXTRA
% -------------------------------------------------------------------------
% [ ] Investigate regularisation and iteration settings
% [ ] Work on nii.gz? (see: https://uk.mathworks.com/matlabcentral/fileexchange/47698-savezip)
%
% -------------------------------------------------------------------------
% VALIDATION
% -------------------------------------------------------------------------
% [ ] https://www.synapse.org/#!Synapse:syn3207203 
% [ ] fmriproc
%
% -------------------------------------------------------------------------
% POPULATIONS
% -------------------------------------------------------------------------
%     | Name         | Modality     | Labels                                   | NumSubj 
% -------------------------------------------------------------------------
% 1   | ATLAS        | T1           | 1.les                                    | 142     
% 2   | BALGRIST     | T1,PD        | 1.spn                                    | 19      
% 3   | CROMIS       | CT           | 1.bg,2.bone                              | 686     
% 4   | CROMISLABELS | CT           | 1.les,2.cal                              | 60      
% 5   | DELIRIUM     | CT           | n/a                                      | 1,025   
% 6   | IXI          | T1,T2,PD,MRA | n/a                                      | 567     
% 7   | IXIRC        | GM,WM,CSF    | n/a                                      | 32      
% 8   | IXIC         | GM,WM,CSF    | n/a                                      | 32      
% 9   | MICCAI2012   | T1           | 1.cgm,2.sgm,3.spn,4.wm,5,csf,6.ven       | 35     
% 10  | MRBRAINS18   | T1           | 1.cgm,2.sgm,3.spn,4.wm,5.cer,6.csf,7.ven | 7       
% 11  | ROB          | CT           | n/a                                      | 72    
% 12  | MPMCOMPLIANT | MPM          | n/a                                      | 10
% 13  | MADRID       | T1           | n/a                                      | 16
% 14  | IBSR18       | T1           | n/a                                      | 18
% 15  | CTHEALTHY    | CT           | 1.bg,2.bone                              | 50
% 16  | LPBA40       | T1           | n/a                                      | 40
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

%%%%%%%%%%%%%%%%%%%
% Parse input struct
%------------------

if nargin < 1, opt = struct; end

% For setting paths and dirs
if ~isfield(opt,'user'), opt.user = 'mbrud-home'; end 
% Fit on 3D or 2D data?
if ~isfield(opt,'run3d'), opt.run3d = false; end 
% 2D plane ['ax','cor','sag']
if ~isfield(opt,'ax2d'), opt.ax2d = 'ax'; end 
% Number of subjects
if ~isfield(opt,'numsubj'), opt.numsubj = Inf; end           
% Model 0 | Testing
% Model 1 | Labels are used (K1=10), trying to get nice GM, WM and CSF
% Model 2 | Labels are not used (K1=12), unsupervised for better normalisation
% Model 3 | Labels are used (K1=7, mg_ix=2), trying to get nice GM, WM and CSF
% Model 4 | Fit T1 (use to init) and CT (K1=12), unsupervised
% Model 5 | Fit a learned model to new subjects
if ~isfield(opt,'models'), opt.models = 0; end     
% Number of parfor workers
if ~isfield(opt,'nw'),     opt.nw = 0; end   
% Show figures
if ~isfield(opt,'show'),   opt.show = true; end   

%%%%%%%%%%%%%%%%%%%
% Set/get user specific
%------------------

[dir_data,dir_res] = UserSpecific(opt.user);

%%%%%%%%%%%%%%%%%%%
% Define populations
%------------------
  
[P,ix] = GetPopulations(dir_data,opt);

for model_num=opt.models % loop over models
        
    %%%%%%%%%%%%%%%%%%%
    % Get Model
    %------------------

    [P1,sett,model] = GetModel(model_num,P,ix,dir_res,opt);

    if ~isempty(P1)
        %%%%%%%%%%%%%%%%%%%
        % Fit model
        %------------------

        if isempty(model)
            % Groupwise registration
            FitModel('groupwise',P1,sett); 
        else
            % Fit already learned model
            FitModel('register',P1,model,sett); % 'model' is loaded in function GetModel()

        %     % Evaluate model
        %     res = GetGoodnessOfFit(sett.write.dir_res);
        end

        if 0
           spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(wimc).*\.nii$')); 
        %    spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(wc|wimc).*\.nii$')); 
        %    spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(im[1,2,3]).*\.nii$')); 
        end
    end
end
end
%==========================================================================

%==========================================================================
function [P,ix] = GetPopulations(dir_data,opt)
% Get available populations
%
% P{i} = {'NAME',{MOD1,..,MODC},NumSubj,PopIx,cm_map,CT,do_bf,do_dc}

ax2d  = opt.ax2d;
run3d = opt.run3d;

if run3d
    % 3D data will be used
    d_2D = '';
else
    % Fit on 2D data, anatomical plane decided by AX2D variable
    fprintf('=============\nOBS! USING 2D DATA!\n=============\n\n');    
    d_2D = fullfile('2D',ax2d);
end

ix = struct('ATLAS',1,'BALGRIST',2,'CROMIS',3,'CROMISLABELS',4,'DELIRIUM',5, ...
            'IXI',6,'IXIC',7,'IXIRC',8,'MICCAI2012',9,'MRBRAINS18',10,'ROB',11, ...
            'MPMCOMPLIANT',12,'MADRID',13,'IBSR18',14,'CTHEALTHY',15, ...
            'LPBA40',16);
P  = cell(1,numel(ix));

% Populations of images (GMM will be fitted)
P{ix.ATLAS}        = {fullfile(dir_data,d_2D,'ATLAS'),        {'T1'}, Inf, [], {}, false, true, true};
P{ix.BALGRIST}     = {fullfile(dir_data,d_2D,'BALGRIST'),     {'T1','PD'}, Inf, [], {}, false, true, true};
P{ix.CROMIS}       = {fullfile(dir_data,d_2D,'CROMIS'),       {'CT'}, Inf, [], {}, true, true, false};
P{ix.CROMISLABELS} = {fullfile(dir_data,d_2D,'CROMISLABELS'), {'CT'}, Inf, [], {}, true, true, false};
P{ix.DELIRIUM}     = {fullfile(dir_data,d_2D,'DELIRIUM'),     {'CT'}, Inf, [], {}, true, true, false};
P{ix.IXI}          = {fullfile(dir_data,d_2D,'IXI'),          {'T1','T2','PD','MRA'}, Inf, [], {}, false, true, true};
P{ix.MICCAI2012}   = {fullfile(dir_data,d_2D,'MICCAI2012'),   {'T1'}, Inf, [], {}, false, true, true};
P{ix.MRBRAINS18}   = {fullfile(dir_data,d_2D,'MRBRAINS18'),   {'T1'}, Inf, [], {}, false, true, true};
P{ix.ROB}          = {fullfile(dir_data,d_2D,'ROB'),          {'CT'}, Inf, [], {}, true, true, false};
P{ix.MPMCOMPLIANT} = {fullfile(dir_data,d_2D,'MPMCOMPLIANT'), {'MT','PD','R2','T1'}, Inf, [], {}, false, true, [true true false true]};
P{ix.MADRID}       = {fullfile(dir_data,d_2D,'MADRID'),       {'T1'}, Inf, [], {}, false, true, true};
P{ix.IBSR18}       = {fullfile(dir_data,d_2D,'IBSR18'),       {'T1'}, Inf, [], {}, false, true, true};
P{ix.CTHEALTHY}    = {fullfile(dir_data,d_2D,'CTHEALTHY'),    {'CT'}, Inf, [], {}, true, true, false};
P{ix.LPBA40}       = {fullfile(dir_data,d_2D,'LPBA40'),       {'T1'}, Inf, [], {}, false, true, true};

% Populations of tissue segmentations (2D not available)
P{ix.IXIC}  = {fullfile(dir_data,'IXIC'),  {'GM','WM','CSF'}, Inf, [], {}, false};
P{ix.IXIRC} = {fullfile(dir_data,'IXIRC'), {'GM','WM'}, Inf, [], {}, false};    
end
%==========================================================================

%==========================================================================
function [dir_data,dir_res] = UserSpecific(user)
% Sets MATLAB path to diffeo-segment and auxiliary-functions, and return
% path to data directory and where to write output results
%
% Available users
% . mbrud-home
% . mbrud-fil
if strcmp(user,'mbrud-home')    
    addpath('/home/smajjk/dev/diffeo-segment')
    addpath('/home/smajjk/dev/auxiliary-functions')
    dir_data = '/home/smajjk/Data/Nii/diffeo-segment/';
    dir_res  = '/home/smajjk/Data/Results/diffeo-segment';
elseif strcmp(user,'mbrud-fil')    
    addpath('/home/mbrud/dev/mbrud/code/matlab/diffeo-segment')
    addpath('/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions')
    dir_data = '/scratch/Nii/TrainingData/diffeo-segment/';
    dir_res  = '/scratch/Results/diffeo-segment';
elseif strcmp(user,'yael-fil')    
    addpath('~/Devel/balbasty/auxiliary-functions/');
    addpath('~/Devel/computational-anatomy/diffeo-segment-pca/');
    dir_data = '/scratch/new_segment/TrainingData';
    dir_res  = '/export/data/experiments/new_segment/2020/';
elseif strcmp(user,'yael-fil-master')    
    addpath('~/Devel/balbasty/auxiliary-functions/');
    addpath('~/Devel/computational-anatomy/diffeo-segment/');
    dir_data = '/scratch/new_segment/TrainingData';
    dir_res  = '/export/data/experiments/new_segment/200110/';
elseif strcmp(user,'gold-fil')    
    addpath('/scratch/gold/mikael/code/auxiliary-functions/');
    addpath('/scratch/gold/mikael/code/diffeo-segment/');
    dir_data = '/scratch/gold/mikael/data';
    dir_res  = '/scratch/gold/mikael/results';    
else
    error('Undefined user!');
end
end
%==========================================================================