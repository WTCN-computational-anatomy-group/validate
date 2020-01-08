function FindBestModel
% Find model with best fit

dir_exp     = 'experiment';
experiments = 0:6;

mx_mce = -Inf;
mx_mds = -Inf;

for e=experiments
    d = fullfile(dir_exp,num2str(e));
    f = spm_select('FPList',d,'^res.*\.mat$');
    
    fi = deblank(f(e + 1,:));
    load(fi);
        
%     mce = median(res.ce);
%     mds = median(res.ds);
%     mds
%     if mds > mx_mds
%         mx_mds = mds;
%         emx    = e;
%     end
    
    for i=1:size(f,1)
        fi = deblank(f(i,:));
        load(fi);
        
        mce = median(res.ce);
        mds = median(res.ds);
        
        if mds > mx_mds
            mx_mds = mds;
            emx    = e;
            imx    = i;
            tmp = sort(res.ds);
        end
    end
end
end
%==========================================================================