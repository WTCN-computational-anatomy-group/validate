function res = GetGoodnessOfFit(dir_lab)
% Evaluate model fit

% Get normalised labels
f   = spm_select('FPList',dir_lab);
Nii = nifti(f);
N   = numel(Nii);

% Compute dice score
cl  = zeros(1,0.5*N*(N - 1));
res = struct('ce',cl,'ds',cl);
cnt = 1;
for i1=1:N
    im1 = Nii(i1).dat(:);
    for i2=i1+1:numel(Nii)
        im2 = Nii(i2).dat(:);
                
        % Compute average cross entropy
        ce = -sum(im1.*log(im2 + eps)) - sum(im2.*log(im1 + eps));
        
        % Compute dice score
        ds = dice(im1,im2);
        
        res.ce(cnt) = ce;
        res.ds(cnt) = ds;        
        cnt         = cnt + 1;
    end
end

if false
    figure(666)
    subplot(121)
    boxplot(res.ce)
    grid on
    axis tight
    title('CE')
    subplot(122)
    boxplot(res.ds)
    grid on
    axis tight
    title('DS')
end
end
%==========================================================================