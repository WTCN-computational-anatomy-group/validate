clear;

DO_BRAINWEB = 2; % 0: not, 1: 2D, 2: 3D
DO_RIRE     = false;

if DO_BRAINWEB
    % Simulated data
    %-----------------------

    DirData = 'test-data/BrainWeb';

    for i=1:10
        % Simulate some misaligned BrainWeb images (T1, T2, PD)        
        [Nii3d,Nii2d] = Simulate('DirRef',DirData, ...
                                 'Random',[true true false, ... % do_sim, add_bias, misalign
                                          false false false]);   % thick_slice, crop_fov, add_noise
    
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
    
    SHOW_ALIGNED = true;
    DO_PW        = false;
    
    % Get images
    DirData      = '/home/mbrud/Data/Challenges/RIRE/train/training_001';
%     DirData      = '/home/mbrud/Data/Challenges/RIRE/test/patient_001';
    [Nii,info,C] = LoadRIRE(DirData);    
    IxFixed      = info.t1.ix;
    scans        = [info.t1.ix info.t2.ix info.pd.ix info.ct.ix info.pet.ix];
    
    if DO_PW
        % Select scans of interest
        IxFixed = info.t1.ix;
        scans   = [IxFixed info.ct.ix];              
    end
    
    scans = sort(scans);
    Nii   = Nii(scans);
        
    % Registration options
    IxFixed = find(scans == IxFixed);
    opt     = struct('IxFixed', IxFixed, ...
                     'ShowAlign', 1, ...
                     'Samp', [8 4 2 1], ...
                     'PreComp', true, ...
                     'ShowFit4Scaling', false, ...
                     'DegBoundCond', [2 0]);
    
    % Make results dir    
    DirRes = fullfile(DirData,'results');
    if ~isfolder(DirRes), mkdir(DirRes); end
    
    % Do registration             
    tic
    [~,R] = spm_njtv_coreg(Nii,opt);
    toc
        
    % Make temp folder
    DirCpy = 'temp';
    if isfolder(DirCpy), rmdir(DirCpy,'s'); end; mkdir(DirCpy);

    % Make C transformations in R    
    R0            = repmat(eye(4),[1 1 C]);
    R0(:,:,scans) = R;
    
    % Save registration results
    save(fullfile(DirRes,'R.mat'),'R0');
                  
    % Apply transformation (to copies) and save in 'DirCpy'
    nf = cell(1,numel(Nii));
    for c=1:numel(Nii)
       Mmov        = Nii(c).mat;        
       f           = Nii(c).dat.fname;
       [~,nam,ext] = fileparts(f);
       nf{c}       = fullfile(DirCpy,[nam ext]);
       copyfile(f, nf{c})

       M = R(:,:,c)*Mmov;
       spm_get_space(nf{c},M);
    end

    if contains(DirData, 'train')
        % Get error
        ValidateRIRE(DirData,DirRes);
    end
    
    if SHOW_ALIGNED
        % Show aligned copies
        spm_check_registration(char(nf));    
    end
end