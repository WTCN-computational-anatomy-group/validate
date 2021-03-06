function MakeTrainingData
%
% -------------------------------------------------------------------------
% POPULATIONS
% -------------------------------------------------------------------------
% ATLAS         | T1           | yes | N = 142
% BALGRIST      | T1           | yes | N = 19
% CROMIS        | CT           | n/a | N = 626
% CROMISLABELS  | CT           | yes | N = 60
% CROMISPETTERI | CT           | yes | N = 20
% DELIRIUM      | CT           | n/a | N = 1,025
% IXI           | T1,T2,PD,MRA | n/a | N = 567
% MADRID        | T1           | n/a | N = 16
% MICCAI2012    | T1           | yes | N = 35
% MPMCOMPLIANT  | MPM          | n/a | N = 10
% MRBRAINS18    | T1           | yes | N = 7 
% RIRE          | T1,T2,PD,CT  | n/a | N = 19
% ROB           | CT           | n/a | N = 72
% -------------------------------------------------------------------------
%
%__________________________________________________________________________


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global settings
%-------------------------------------

S0         = 8;
NumWorkers = Inf;
DirData0   = '/media/mbrud/Storage/Original';
DirOut     = '/home/mbrud/Data/Training/diffeo-segment-new';
Write2D    = false;

if ~(exist(DirOut,'dir') == 7), mkdir(DirOut); end

Do = false;
Do = struct('ATLAS',Do,'BALGRIST',Do,'CROMIS',Do,'CROMISLABELS',Do, ...
            'CROMISPETTERI',Do,'DELIRIUM',Do,'IXI',Do,'MADRID',Do,'MICCAI2012',Do, ...
            'MPMCOMPLIANT',Do,'MRBRAINS18',Do,'RIRE',Do, 'ROB',Do);   
        
% Do.ATLAS         = true;
% Do.BALGRIST      = true;
% Do.CROMIS        = true;
% Do.CROMISLABELS  = true;
% Do.CROMISPETTERI = true;
% Do.DELIRIUM      = true;
% Do.IXI           = true;
% Do.MADRID        = true;
Do.MICCAI2012    = true;
% Do.MPMCOMPLIANT  = true;
% Do.MRBRAINS18    = true;
% Do.RIRE          = true;
% Do.ROB           = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATLAS
%-------------------------------------
Population = 'ATLAS';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:2:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:))); % T1
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 1,:))); % Labels

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni = true; 
    opt.do.crop     = true;
    
    opt.do.vx          = true;
    opt.crop.keep_neck = false;
    
    opt.do.nm_reorient = true;
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALGRIST
%-------------------------------------
Population = 'BALGRIST';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    % Get files
    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    % Build dat
    dat = struct; cnt = 1;
    for s=1:3:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));     % T1    
        dat(cnt).Nii{1}(2) = nifti(deblank(files(s + 1,:))); % PD    
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 2,:))); % Labels
        dat(cnt).Nii{2}(2) = nifti;

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni = true;
    opt.do.crop     = true;
    opt.do.coreg    = true;
    opt.do.reslice  = true;
    opt.reslice.ref = 1;
    opt.crop.keep_neck = false;
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROMIS
%-------------------------------------
Population = 'CROMIS';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;    
    opt.dir_out     = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig = true;
    opt.do.real_mni = true;
    opt.do.crop     = true;
    opt.crop.keep_neck = false;
    opt.do.vx       = true; 
    if Write2D
        opt.do.nm_reorient = true;        
        opt.do.write2d     = true;
        opt.vx.min_1mm     = false;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROMISLABELS
%-------------------------------------
Population = 'CROMISLABELS';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:2:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:))); % CT
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 1,:))); % Labels

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;    
    opt.dir_out     = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig = true;
    opt.do.real_mni = true;
    opt.do.crop     = true;
    opt.crop.keep_neck = false;
    opt.do.vx       = true; 
    opt.labels.part = {1,2};
    if Write2D
        opt.do.nm_reorient = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROMISPETTERI
%-------------------------------------
Population = 'CROMISPETTERI';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:2:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:))); % CT
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 1,:))); % Labels

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt                = struct;    
    opt.dir_out        = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig    = true;
    opt.do.real_mni    = true;
    opt.do.crop        = true;    
    opt.crop.keep_neck = false;
%     opt.do.vx          = true;
%     opt.do.nm_reorient = true;
    if Write2D         
        opt.do.vx          = true;
%         opt.do.nm_reorient = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=7, RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELIRIUM
%-------------------------------------
Population = 'DELIRIUM';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;        
    opt.dir_out     = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig = true;
    opt.do.real_mni = true;
    opt.do.crop     = true;
    opt.crop.keep_neck = false;
    opt.do.vx       = true;     
    if Write2D
        opt.do.nm_reorient = true;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IXI
%-------------------------------------
Population = 'IXI';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:4:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));     % MRA
        dat(cnt).Nii{1}(2) = nifti(deblank(files(s + 1,:))); % PD
        dat(cnt).Nii{1}(3) = nifti(deblank(files(s + 2,:))); % T1
        dat(cnt).Nii{1}(4) = nifti(deblank(files(s + 3,:))); % T2

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni    = true; 
    opt.do.coreg       = true;
    opt.do.crop        = true;
    opt.crop.keep_neck = false; 
    opt.do.reslice     = true;
    opt.reslice.ref    = 3;
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;        
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MADRID
%-------------------------------------
Population = 'MADRID';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:))); % T1

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni = true; 
    opt.do.crop     = true;
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MICCAI2012
%-------------------------------------
Population = 'MICCAI2012';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:2:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:))); % T1
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 1,:))); % Labels

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);   
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end    
    opt.do.real_mni = true;    
    opt.do.crop     = true;        
    opt.do.erode    = true;
    opt.crop.keep_neck = false;
    opt.labels.part = {[23 30 31 32 36 37 47 48 55 56 57 58 59 60 61 62 75 76 38 39 71 72 73 100 101 102 103 104 105 106 107 108 109 112 113 114 115 116 117 118 119 120 121 122 123 124 125 128 129 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207], ...                    
                       [35 40 41 44 45 69]};
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers <= 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPMCOMPLIANT
%-------------------------------------
Population = 'MPMCOMPLIANT';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:4:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));     % MT
        dat(cnt).Nii{1}(2) = nifti(deblank(files(s + 1,:))); % PD
        dat(cnt).Nii{1}(3) = nifti(deblank(files(s + 2,:))); % R2
        dat(cnt).Nii{1}(4) = nifti(deblank(files(s + 3,:))); % T1

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni = true;
    opt.do.coreg    = true;
    opt.do.crop     = true;
    opt.do.reslice  = true;
    opt.reslice.ref = 4;
    if Write2D        
        opt.crop.keep_neck = false;
        opt.do.nm_reorient = true;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRBRAINS18
%-------------------------------------
Population = 'MRBRAINS18';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:4:S
        dat(cnt).Nii{1}    = nifti;
%         dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));     % FLAIR
%         dat(cnt).Nii{1}(2) = nifti(deblank(files(s + 1,:))); % IR
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s + 2,:))); % T1
        dat(cnt).Nii{2}    = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 3,:))); % Labels
%         dat(cnt).Nii{2}(2) = nifti;
%         dat(cnt).Nii{2}(3) = nifti;

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt                = struct;
    opt.dir_out        = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni    = true; 
    opt.do.crop        = true;
    opt.crop.keep_neck = false;
%     opt.do.coreg       = true;
%     opt.do.reslice     = true;    
%     opt.reslice.ref    = 3;
    opt.labels.part    = {1,2,8,3,7,5,6}; % cgm,sgm,spn,wm,cer,csf,ven
    opt.do.skullstrip  = true;
    opt.do.bfcorr      = true;
    opt.do.bb_spm      = true;
    if Write2D
        opt.do.nm_reorient = true;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIRE
%-------------------------------------
Population = 'RIRE';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    [~,d0]   = spm_select('FPList',DirData,'^.*\.nii$');
    S       = size(d0,1);

    dat = struct; cnt = 1;
    for s=1:S
        
        d        = deblank(d0(s,:));
        files_ct = spm_select('FPListRec',d,'_ct.nii$');
        files_t1 = spm_select('FPListRec',d,'_mr_T1.nii$');
        files_t2 = spm_select('FPListRec',d,'_mr_T2.nii$');
        files_pd = spm_select('FPListRec',d,'_mr_PD.nii$');
    
        if isempty(files_ct) || isempty(files_t1) || isempty(files_t2) || isempty(files_pd)
            continue
        end
            
        dat(cnt).Nii{1}    = nifti;           
        dat(cnt).Nii{1}(1) = nifti(deblank(files_t1));
        dat(cnt).Nii{1}(2) = nifti(deblank(files_t2));
        dat(cnt).Nii{1}(3) = nifti(deblank(files_pd));
        dat(cnt).Nii{1}(4) = nifti(deblank(files_ct));
        
        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;
    opt.dir_out     = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni = true; 
    opt.do.crop     = true;
    opt.do.coreg    = true;
    opt.do.reslice  = true;
    opt.reslice.ref = 1;
    opt.do.vx       = true;
    if Write2D
        opt.do.nm_reorient = true;
        opt.crop.keep_neck = false;
        opt.do.vx          = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROB
%-------------------------------------
Population = 'ROB';
if Do.(Population)
    fprintf('======================================\n')
    fprintf('| %s\n',Population);
    fprintf('======================================\n')

    DirData = fullfile(DirData0,Population);
    files   = spm_select('FPListRec',DirData,'^.*\.nii$');
    S       = size(files,1);

    dat = struct; cnt = 1;
    for s=1:S
        dat(cnt).Nii{1}    = nifti;
        dat(cnt).Nii{1}(1) = nifti(deblank(files(s,:)));

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt             = struct;        
    opt.dir_out     = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end    
    opt.do.real_mni = true;
    opt.do.crop     = true;
    opt.crop.keep_neck = false;
    opt.do.vx       = true; 
    if Write2D
        opt.do.nm_reorient = true;
        opt.vx.min_1mm     = false;
        opt.do.write2d     = true;
        opt.dir_out2d      = fullfile(DirOut,['2D_' Population]);
        if (exist(opt.dir_out2d,'dir') == 7), rmdir(opt.dir_out2d,'s'); end
    end
    opt.do.go2native = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), RunPreproc(dat(s).Nii,opt); end
    end
end