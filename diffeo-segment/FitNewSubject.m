clear;

% Add required toolboxes to MATLAB's path
addpath('/home/mbrud/dev/mbrud/code/matlab/diffeo-segment')      % https://github.com/WTCN-computational-anatomy-group/diffeo-segment
addpath('/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions') % https://github.com/WTCN-computational-anatomy-group/auxiliary-functions

DirModel = '/scratch/Results/diffeo-segment/results/model-3-2D-ax';
DirData  = '/scratch/Nii/TrainingData/diffeo-segment/2D/ax/MICCAI2012';
DirOut   = '/scratch/Results/diffeo-segment/results/FitNewSubject';
PthModel = fullfile(DirModel,'model_spm_mb.mat');

% Data
pths_im  = spm_select('FPList',DirData,'^((?!_glm).)*\.nii$');
pths_lab = spm_select('FPList',DirData,'_glm.*\.nii$');
            
N    = size(pths_im,1);
dat0 = struct('F',[],'labels',[]);
for n=1:N
    im = nifti(deblank(pths_im(n,:)));
    lab = {nifti(deblank(pths_lab(n,:))),{}};
    if n == 1
        dat        = dat0;
        dat.F      = im;
        dat.labels = lab;
    else
        datn        = dat0;
        datn.F      = im;
        datn.labels = lab;
        dat         = [dat datn];
    end
end

% Pick indices
n   = 21;
dat = dat(n);

% Settings
sett                   = struct;
sett.show.figs         = {'model'};%,'segmentations','normalised','parameters','intensity'};
sett.write.dir_res     = DirOut;
sett.write.tc          = [true false false];
sett.write.im          = [true false false false];
sett.write.df          = [false false];
sett.write.labels      = [false false];
sett.do.infer          = true;
sett.show.mx_subjects  = 4;
% sett.gen.samp_min      = 3;
% sett.clean_z.mrf       = 2;
% sett.clean_z.gwc_tix   = struct('gm',[1],'wm',[2],'csf',[3]);
sett.write.vel         = true;
sett.write.affine      = true;

% Run Register
[dat,mu,sett] = spm_mb_fit(dat,'PthModel',PthModel,'sett',sett);

% Write results in normalised space
spm_mb_output(dat,mu,sett);

spm_check_registration(spm_select('FPList',sett.write.dir_res,'^(im1|c).*\.nii$'))