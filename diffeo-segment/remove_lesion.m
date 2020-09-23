clear

N = 25;

% ATLAS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/ATLAS';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);


Nii_im = nifti(pths_im);
Nii_lab = nifti(pths_lab);
for n=1:N
    im = Nii_im(n).dat();
    mask = Nii_lab(n).dat() == 1;
    im(mask) = 0;
    Nii_im(n).dat(:,:,:) = im;
    
    lab = Nii_lab(n).dat();
    lab(mask) = 0;
    Nii_lab(n).dat(:,:,:) = lab;
end

%%

N = 25;

% CROMIS
dir_data = '/home/mbrud/Data/Training/diffeo-segment/CROMISPETTERI';
pths_im  = spm_select('FPList',dir_data,'^((?!_labels).)*\.nii$');
pths_lab = spm_select('FPList',dir_data,'_labels.*\.nii$');        
pths_im  = pths_im(1:min(N,size(pths_im, 1)), :);
pths_lab = pths_lab(1:min(N,size(pths_lab, 1)), :);    
pths_im  = cellstr(pths_im);
pths_lab = cellstr(pths_lab);


Nii_im = nifti(pths_im);
Nii_lab = nifti(pths_lab);
for n=1:N
    
    lab = Nii_lab(n).dat();
    lab(lab == 5) = 4;
    lab(lab == 6) = 5;
    lab(lab == 7) = 6;
    Nii_lab(n).dat(:,:,:) = lab;
end