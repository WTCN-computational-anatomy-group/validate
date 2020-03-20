clear;

DO_BRAINWEB = 1; % 0: not, 1: 2D, 2: 3D
DO_RIRE     = true;

if DO_BRAINWEB
    % Simulated data
    %-----------------------

    DirData = 'test-data/BrainWeb';

    for i=1:10
        % Simulate some misaligned BrainWeb images (T1, T2, PD)        
        [Nii3d,Nii2d] = Simulate('DirRef',DirData, ...
                                 'Random',[true true true, ... % do_sim, add_bias, misalign
                                          true true true]);   % thick_slice, crop_fov, add_noise
    
        % Do registration             
        tic
        if     DO_BRAINWEB == 1, spm_njtv_coreg(Nii2d); % Run on 2D images
        elseif DO_BRAINWEB == 2, spm_njtv_coreg(Nii3d); % Run on 3D images
        end
        toc
    end
end

if DO_RIRE    
    % Hospital data
    %-----------------------
    
    DirData = 'test-data/RIRE-training001';
    Nii     = nifti(spm_select('FPList',DirData,'^.*\.nii$'));    
    Nii     = Nii([3 1 2 4 5]); % t1 ct pd t2 pet 3 1 2 4 5
    
    % Do registration             
    tic
    [~,R,ix] = spm_njtv_coreg(Nii);
    toc
    
    ixf = ix.fixed;
    ixm = ix.moving;
    
    % Apply transformation (to copies) and display
    DirCpy = 'temp';
    if isfolder(DirCpy), rmdir(DirCpy,'s'); end; mkdir(DirCpy);
            
    save(fullfile(DirCpy,'R.mat'),'R','ixf','ixm');
    
    nf = cell(1,numel(Nii));
    for c=1:numel(Nii)
       f           = Nii(c).dat.fname;
       [~,nam,ext] = fileparts(f);
       nf{c}       = fullfile(DirCpy,[nam ext]);
       copyfile(f, nf{c})
       
       if c == ixf, continue; end
       
       Mmov = Nii(c).mat;
       M    = R(:,:,c)\Mmov;
       spm_get_space(nf{c},M);
    end
    
    spm_check_registration(char(nf))
end
