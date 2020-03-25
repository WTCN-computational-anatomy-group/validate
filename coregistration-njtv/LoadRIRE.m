function [Nii,info,C] = LoadRIRE(DirData)
Files = spm_select('FPListRec',DirData,'^.*\.nii$');
Nii   = nifti(Files);
C     = size(Files,1);
s     = struct('ix',0,'nam','');
info  = struct('t1',s,'t1r',s,'t2',s,'t2r',s,'pd',s,'pdr',s,'ct',s,'pet',s,'mpr',s);
for c=1:C
    fn = deblank(Files(c,:));
    d  = strsplit(fn,filesep);
    d  = d{end - 1};
    
    if strcmp(d,'mr_T1')
        info.t1.ix  = c; 
        info.t1.nam = 'T1';
    end
    if strcmp(d,'mr_T1_rectified')
        info.t1r.ix  = c;
        info.t1r.nam = 'T1_rectified';
    end
    if strcmp(d,'mr_T2')
        info.t2.ix  = c; 
        info.t2.nam = 'T2';
    end
    if strcmp(d,'mr_T2_rectified')
        info.t2r.ix  = c; 
        info.t2r.nam = 'T2_rectified';
    end
    if strcmp(d,'mr_PD')
        info.pd.ix  = c; 
        info.pd.nam = 'PD';
    end
    if strcmp(d,'mr_PD_rectified')
        info.pdr.ix  = c; 
        info.pdr.nam = 'PD_rectified';
    end
    if strcmp(d,'ct')
        info.ct.ix  = c; 
        info.ct.nam = 'ct';
    end
    if strcmp(d,'pet') 
        info.pet.ix  = c; 
        info.pet.nam = 'pet';
    end
    if strcmp(d,'mr_MP-RAGE') 
        info.mpr.ix  = c; 
        info.mpr.nam = 'MP-RAGE';
    end
end

fn = fieldnames(info);
for i=1:numel(fn)
    if info.(fn{i}).ix == 0
        info = rmfield(info,fn{i});
    end
end
C = numel(fieldnames(info));
%==========================================================================