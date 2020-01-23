addpath('/home/mbrud/dev/mbrud/code/matlab/Patient-Preprocessing/');
%
% POPULATIONS
%
% [x] | 1.ATLAS        | T1           | 1.les                 | N = 142
% [x] | 2.BALGRIST     | T1           | 1.spn                 | N = 19
% [x] | 3.CROMIS       | CT           | n/a                   | N = 686
% [x] | 4.CROMISLABELS | CT           | 1.les,2.cal           | N = 60
% [x] | 5.DELIRIUM     | CT           | n/a                   | N = 1,025
% [x] | 6.IXI          | T1,T2,PD,MRA | n/a                   | N = 567
% [x] | 7.MICCAI2012   | T1           | 1.gm,2.wn,3.ven       | N = 35
% [x] | 8.MRBRAINS18   | T1,FLAIR,IR  | 1.gm,2.wn,3.ven,4.cer | N = 7 
% [ ] | 9.RIRE         | T1,T2,PD,CT  | n/a                   | N = 19
% [x] | 10.ROB         | CT           | n/a                   | N = 72
%
%__________________________________________________________________________

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global settings
%-------------------------------------

S0         = Inf;
NumWorkers = 8;
DirData0   = '/scratch/Nii/Original';
DirOut     = '/scratch/Nii/TrainingData/diffeo-segment';

if ~(exist(DirOut,'dir') == 7), mkdir(DirOut); end

Do = false;
Do = struct('ATLAS',Do,'BALGRIST',Do,'CROMIS',Do,'CROMISLABELS',Do, ...
            'DELIRIUM',Do,'IXI',Do,'MICCAI2012',Do,'MRBRAINS18',Do,'RIRE',Do, ...
            'ROB',Do);   
        
% Do.ATLAS      = true;
% Do.BALGRIST   = true;
% Do.DELIRIUM   = true;
% Do.IXI        = true;
% Do.MICCAI2012 = true;
% Do.MRBRAINS18 = true;
% Do.ROB = true;

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
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true; 
    opt.do.crop         = true;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);     
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALGRIST
%-------------------------------------
Population = 'BALGRIST';
if Do.BALGRIST
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
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.do.coreg        = true;
    opt.do.reslice      = true;
    opt.reslice.ref     = 1;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt         = struct;
    opt.do.vx   = true; 
    opt.vx.size = 1;        
    opt.vx.deg  = 0;
    opt.dir_out         = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig     = true;
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt         = struct;
    opt.do.vx   = true; 
    opt.vx.size = 1;        
    opt.vx.deg  = 0;
    opt.dir_out         = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig     = true;
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt         = struct;
    opt.do.vx   = true; 
    opt.vx.size = 1;        
    opt.vx.deg  = 0;
    opt.dir_out         = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.res_orig     = true;
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true; 
    opt.do.coreg        = true;
    opt.do.crop         = true;
    opt.do.reslice      = true;
    opt.reslice.ref     = 3;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);   
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.labels.part     = {[23 30 31 32 36 37 47 48 55 56 57 58 59 60 61 62 75 76 38 39 71 72 73 100 101 102 103 104 105 106 107 108 109 112 113 114 115 116 117 118 119 120 121 122 123 124 125 128 129 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207], ...
                           [23 30 31 32 36 37 47 48 55 56 57 58 59 60 61 62 75 76 38 39 71 72 73], ...
                           [35], ...
                           [40 41 44 45 69], ...
                           46, ...
                           [4 11 49 50 51 52]};
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers <= 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
%         dat(cnt).Nii{2}(1) = nifti;
        dat(cnt).Nii{2}(1) = nifti(deblank(files(s + 3,:))); % Labels
%         dat(cnt).Nii{2}(3) = nifti;

        if cnt == S0, break; end
        cnt = cnt + 1;    
    end

    % Set options and do preprocessing
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true; 
    opt.do.crop         = true;
    opt.do.coreg        = false;
    opt.do.reslice      = false;
    opt.reslice.ref     = 1;
    opt.labels.part     = {1,2,8,3,7,5,6}; % cgm,sgm,spn,wm,cer,csf,ven
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt = struct;
    opt.dir_out         = fullfile(DirOut,Population);     
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end
    opt.do.real_mni     = true; 
    opt.do.crop         = true;
    opt.do.coreg        = true;
    opt.do.reslice      = true;
    opt.reslice.ref     = 3;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
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
    opt         = struct;
    opt.dir_out         = fullfile(DirOut,Population);    
    if (exist(opt.dir_out,'dir') == 7), rmdir(opt.dir_out,'s'); end    
    opt.do.real_mni     = true;
    opt.do.crop         = true;
    opt.do.write2d      = true;
    opt.dir_out2d       = fullfile(DirOut,['2D_' Population]);
    opt.do.go2native    = false;  
    if NumWorkers == 1 || S0 == 1
        for s=1:numel(dat), out = RunPreproc(dat(s).Nii,opt); end
    else
        parfor(s=1:numel(dat),NumWorkers), out = RunPreproc(dat(s).Nii,opt); end
    end
end