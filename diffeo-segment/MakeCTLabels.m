function MakeCTLabels

DirData = {'/scratch/Nii/TrainingData/diffeo-segment/CTHEALTHY'};

bg_mn   = -1020;
bg_mx   = -200;
bone_mn = 100;
bone_mx = 3000;
prfx    = 'lab_';

for i1=1:numel(DirData)
    Nii = nifti(spm_select('FPList',DirData{i1},'^((?!lab_).)*\.nii$'));
    
    N = numel(Nii);
    for i2=1:N
        img  = Nii(i2).dat();
        bone = img >= bone_mn & img <= bone_mx;
        bg   = img >= bg_mn & img <= bg_mx;
        lab  = zeros(size(img),'uint8');
                
        lab(bg)   = 1;
        lab(bone) = 2;
        
        f             = Nii(i2).dat.fname;
        [pth,nam,ext] = fileparts(f);
        nf            = fullfile(pth,[prfx nam ext]);
        
        create_nii(nf,lab,Nii(i2).mat,spm_type('uint8'),'Labels');    
    end
end
%==========================================================================

%==========================================================================
function Nii = create_nii(pth,dat,mat,dtype,descrip,offset,scl_slope,scl_inter)
if nargin<6, offset    = 0; end
if nargin<7, scl_slope = 1; end
if nargin<8, scl_inter = 0; end

if exist(pth,'file')==2, delete(pth); end

Nii         = nifti;
dm          = size(dat);
Nii.dat     = file_array(pth,dm,dtype,offset,scl_slope,scl_inter);
Nii.mat     = mat;
Nii.mat0    = mat;
Nii.descrip = descrip;
create(Nii);

if numel(dm)==4
    Nii.dat(:,:,:,:) = dat;
else
    Nii.dat(:,:,:)   = dat;
end
%==========================================================================