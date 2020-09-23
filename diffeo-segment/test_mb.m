clear

%-----------------------------------------------------------------------
% Job saved on 13-May-2020 16:06:40 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch = {};
matlabbatch{1}.spm.tools.mb.run.mu.create.K = 9;
matlabbatch{1}.spm.tools.mb.run.mu.create.vx = 1.5;
matlabbatch{1}.spm.tools.mb.run.mu.create.mu_settings = [1e-05 0.5 0];
matlabbatch{1}.spm.tools.mb.run.aff = 'SE(3)';
matlabbatch{1}.spm.tools.mb.run.v_settings = [0.0001 0 0.5 0.125 0.5];
matlabbatch{1}.spm.tools.mb.run.onam = 'mb';
matlabbatch{1}.spm.tools.mb.run.odir = {'/home/mbrud/dev/matlab/validate/diffeo-segment/temp/model-w-lesion'};
matlabbatch{1}.spm.tools.mb.run.cat = {{}};
matlabbatch{1}.spm.tools.mb.run.accel = 0.8;
matlabbatch{1}.spm.tools.mb.run.min_dim = 16;
matlabbatch{1}.spm.tools.mb.run.tol = 0.0005;
matlabbatch{1}.spm.tools.mb.run.sampdens = 2;
matlabbatch{1}.spm.tools.mb.run.save = true;
matlabbatch{1}.spm.tools.mb.run.nworker = Inf;

N = 50;
i = 1;

% ATLAS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/ATLAS';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 1;

matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.images = pths_lab;
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.cm_map = {6, [1 2 3 4 5 7 8 9]}; % lesion
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.w = 0.9999;

matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

% CROMIS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/CROMISPETTERI';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');        
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);    
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 2;

matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.images = pths_lab;
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.cm_map = {
                                                             9 % bg
                                                             7 % bone
                                                             4 % ven
                                                             6 % les
                                                             [5 7] % cal
                                                             [1 2 3 4 5] % brain - cal - les - ven
                                                             8 % rest
                                                             [7 8 9]
                                                             }';
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.w = 0.9999;

matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

spm_jobman('run',matlabbatch);

%% -----------------------------------------------------------------------
% Job saved on 13-May-2020 16:06:40 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch = {};
matlabbatch{1}.spm.tools.mb.run.mu.create.K = 8;
matlabbatch{1}.spm.tools.mb.run.mu.create.vx = 1.5;
matlabbatch{1}.spm.tools.mb.run.mu.create.mu_settings = [1e-05 0.5 0];
matlabbatch{1}.spm.tools.mb.run.aff = 'SE(3)';
matlabbatch{1}.spm.tools.mb.run.v_settings = [0.0001 0 0.5 0.125 0.5];
matlabbatch{1}.spm.tools.mb.run.onam = 'mb';
matlabbatch{1}.spm.tools.mb.run.odir = {'/home/mbrud/dev/matlab/validate/diffeo-segment/temp/model-wo-lesion'};
matlabbatch{1}.spm.tools.mb.run.cat = {{}};
matlabbatch{1}.spm.tools.mb.run.accel = 0.8;
matlabbatch{1}.spm.tools.mb.run.min_dim = 16;
matlabbatch{1}.spm.tools.mb.run.tol = 0.0005;
matlabbatch{1}.spm.tools.mb.run.sampdens = 2;
matlabbatch{1}.spm.tools.mb.run.save = true;
matlabbatch{1}.spm.tools.mb.run.nworker = Inf;

N = 25;
i = 1;

% ATLAS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/ATLAS';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 1;
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

% CROMIS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/CROMISPETTERI';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');        
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);    
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 2;

matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.images = pths_lab;
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.cm_map = {
                                                             8 % bg
                                                             6 % bone
                                                             4 % ven
                                                             [5 6] % cal
                                                             [1 2 3 4 5] % brain - cal - les - ven
                                                             7 % rest
                                                             [6 7 8]
                                                             }';
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.w = 0.9999;

matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

spm_jobman('run',matlabbatch);

%% -----------------------------------------------------------------------
% Job saved on 13-May-2020 16:06:40 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch = {};
matlabbatch{1}.spm.tools.mb.run.mu.create.K = 7;
matlabbatch{1}.spm.tools.mb.run.mu.create.vx = 1.5;
matlabbatch{1}.spm.tools.mb.run.mu.create.mu_settings = [1e-05 0.5 0];
matlabbatch{1}.spm.tools.mb.run.aff = 'SE(3)';
matlabbatch{1}.spm.tools.mb.run.v_settings = [0.0001 0 0.5 0.125 0.5];
matlabbatch{1}.spm.tools.mb.run.onam = 'mb';
matlabbatch{1}.spm.tools.mb.run.odir = {'/home/mbrud/dev/matlab/validate/diffeo-segment/temp/model-wo-lesion-K7'};
matlabbatch{1}.spm.tools.mb.run.cat = {{}};
matlabbatch{1}.spm.tools.mb.run.accel = 0.8;
matlabbatch{1}.spm.tools.mb.run.min_dim = 16;
matlabbatch{1}.spm.tools.mb.run.tol = 0.0005;
matlabbatch{1}.spm.tools.mb.run.sampdens = 2;
matlabbatch{1}.spm.tools.mb.run.save = true;
matlabbatch{1}.spm.tools.mb.run.nworker = Inf;

N = 25;
i = 1;

% ATLAS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/ATLAS';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 1;
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

% CROMIS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/CROMISPETTERI';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');        
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);    
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);

matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.images = pths_im;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_reg = 20000;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.inu.inu_co = 40;
matlabbatch{1}.spm.tools.mb.run.gmm(i).chan.modality = 2;

matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.images = pths_lab;
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.cm_map = {
                                                             7 % bg
                                                             5 % bone
                                                             3 % ven
                                                             [1 5] % cal
                                                             [1 2 3 4] % brain - cal - les - ven
                                                             6 % rest
                                                             [5 6 7]
                                                             }';
matlabbatch{1}.spm.tools.mb.run.gmm(i).labels.true.w = 0.9999;

matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.file = {};
matlabbatch{1}.spm.tools.mb.run.gmm(i).pr.update = true;
matlabbatch{1}.spm.tools.mb.run.gmm(i).tol_gmm = 0.0002;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm_miss = 32;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_gmm = 8;
matlabbatch{1}.spm.tools.mb.run.gmm(i).nit_appear = 4;
i = i + 1;

spm_jobman('run',matlabbatch);