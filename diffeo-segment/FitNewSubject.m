clear;

% Add required toolboxes to MATLAB's path
addpath('/home/mbrud/dev/mbrud/code/matlab/diffeo-segment')      % https://github.com/WTCN-computational-anatomy-group/diffeo-segment
addpath('/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions') % https://github.com/WTCN-computational-anatomy-group/auxiliary-functions

% Set paths
DirModel = '/scratch/Results/diffeo-segment/Models/spine-k13-n493';
DirData  = '/scratch/Nii/Original/BALGRIST';
DirOut0  = '/scratch/Results/diffeo-segment/results/FitNewSubject';
PthModel = fullfile(DirModel,'spm_mb_model.mat');

% Get data
Files = spm_select('FPList',DirData,'^((?!spine_labels|_PDw).)*\.nii$');
N     = size(Files,1);
cl    = cell(1,N);
dat   = struct('F',cl);
for n=1:N    
    dat(n).F = nifti(deblank(Files(n,:)));
end

% Pick indices
n    = 13;
datn = dat(n);

% Settings
sett = struct;

sett.show.figs = {'segmentations'};
% sett.show.figs = {'model','segmentations','parameters','intensity'};

DirOut = fullfile(DirOut0,['-n' num2str(n)]);
if isfolder(DirOut), rmdir(DirOut,'s'); end; mkdir(DirOut);

sett.write.dir_res     = DirOut;
sett.write.im          = [true true false false];
sett.write.tc          = false(8,4);
sett.write.tc(1:3,1:3) = true;
sett.write.clean_def   = false;

sett.do.infer = false;

sett.gen.has_spine     = true;
sett.clean_z.mrf       = 2;
sett.clean_z.gwc_tix   = struct('gm',[1],'wm',[2],'csf',[3]);
sett.clean_z.gwc_level = 1;

sett.model.appear_ix  = 1;
sett.model.appear_chn = 1;

sett.nit.init         = 128;
sett.model.init_mu_dm = 8;
sett.var.v_settings   = [0 0 0.2 0.05 0.2]*2;

% Run Register
[datn,mu,sett] = spm_mb_fit(datn,'PthModel',PthModel,'sett',sett);

% Write results in normalised space
spm_mb_output(datn,mu,sett);

spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(imc1|c[1,2,3]).*\.nii$'))
% spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(imc1|c[1,2,3]|wc[1,2,3]).*\.nii$'))