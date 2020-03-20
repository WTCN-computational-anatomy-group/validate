function [oNii,oNii2d,Offset,Rotation,bb0,R3d,R2d] = Simulate(varargin)
% Simulate data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

p              = inputParser;
p.FunctionName = 'SimulateData';
p.addParameter('DirRef','');
p.addParameter('Random',false);
p.addParameter('DirSim','./simulated-data');
p.addParameter('NoisePrct',0);
p.addParameter('Offset',{});
p.addParameter('Rotation',{});
p.addParameter('DownSampling',{});
p.addParameter('DownSamplingDeg',0);
p.addParameter('BiasFieldScl',0);
p.addParameter('ShowSimulated',false);
p.addParameter('ThickSliced',false);
p.addParameter('Margin',{});
p.addParameter('DirPrint','');
p.parse(varargin{:});
DirRef          = p.Results.DirRef;
DirSim          = p.Results.DirSim;
NoisePrct       = p.Results.NoisePrct;     % Percentage of amount of Gaussian noise to add
Offset          = p.Results.Offset;
Rotation        = p.Results.Rotation;
BiasFieldScl    = p.Results.BiasFieldScl;  % Simulate a bias field when value is greater than zero
ShowSimulated   = p.Results.ShowSimulated; % Show simulated images vs reference
mrg             = p.Results.Margin;
DownSampling    = p.Results.DownSampling;
DownSamplingDeg = p.Results.DownSamplingDeg;
DirPrint        = p.Results.DirPrint;
Random          = p.Results.Random;

% Random(1): Do random simulation
% Random(2): In random simulation, include bias field
% Random(3): In random simulation, include rigid transformation
% Random(4): In random simulation, include thick-slices
% Random(5): In random simulation, include smaller FOV
% Random(6): In random simulation, include Rician noise

if 0
    rng('default');
    rand(10);
end

if Random(1)
    %----------------------------------------------------------------------
    % Generate some random parameters
    %----------------------------------------------------------------------
    
    NumChan = 3;
    Idx     = randperm(NumChan);
    Ref     = Idx(1);
    Src     = Idx(Idx ~= Ref);

    % Bias field
    bf_mn = 0;
    bf_mx = 0;
    if numel(Random) >= 2 && Random(2)
        bf_mx = 2;
    end

    t_mn = 0;
    r_mn = 0;
    t_mx = 0;
    r_mx = 0;
    if numel(Random) >= 3 && Random(3)
        % Rotation (degrees)        
        r_mx = 10;
        r_mn = r_mn*pi/180;
        r_mx = r_mx*pi/180;

        % Translation (mm)    
        t_mx = 25;
    end
        
    % Downsampling
    d_mx = 0;
    if numel(Random) >= 4 && Random(4)
        d_mx = 3;
    end    

    % Margin
    Margin = 0;
    if numel(Random) >= 5 && Random(5)
        Margin = 10;
    end
    
    % Bias field
    BiasFieldScl = 0;
    if numel(Random) >= 2 && Random(2)     
        BiasFieldScl = bf_mn + rand*(bf_mx - bf_mn);        
    end

    % Noise
    n_mn = 0;
    n_mx = 0;
    if numel(Random) >= 6 && Random(6)
        n_mx = 0.5;
    end
    NoisePrct = n_mn + rand*(n_mx - n_mn);
        
    Offset   = {[0; 0; 0],[0; 0; 0],[0; 0; 0]};
    Rotation = {[0; 0; 0],[0; 0; 0],[0; 0; 0]};
    for s = Src
        % Translation
        t         = t_mn + (sign(randn(3,1)).*rand(3,1))*(t_mx - t_mn);  
        Offset{s} = t;

        % Rotation
        r           = r_mn + (sign(randn(3,1)).*rand(3,1))*(r_mx - r_mn);
        Rotation{s} = r;
    end

    % Down-sampling
    ds                      = randi(max(d_mx,1),1);
    DownSampling            = {[1 1 1],[1 1 1],[1 1 1]};
    DownSampling{1}(Idx(1)) = 1/ds;
    DownSampling{2}(Idx(2)) = 1/ds;
    DownSampling{3}(Idx(3)) = 1/ds;    
    DownSampling{1}(end)    = 1;
    DownSampling{2}(end)    = 1;
    DownSampling{3}(end)    = 1;
    
    % Margin
    mrg = {[0 0 0],[0 0 0],[0 0 0]};
    if Margin > 0
        mrg1                  = 2*(randperm(NumChan) - 1);
        mrg{1}(Idx(1))        = Margin;
        mrg{1}(1:3 ~= Idx(1)) = mrg1(1);
        mrg{2}(Idx(2))        = Margin;
        mrg{2}(1:3 ~= Idx(2)) = mrg1(2);
        mrg{3}(Idx(3))        = Margin;
        mrg{3}(1:3 ~= Idx(3)) = mrg1(3);
    end
end

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% If true, reslices images to have the same size and orientation matrix.
% Also crops a little bit of the FOV.
Reslice = false;
Padding = -10;

% Determines what plane to extract when creating 2D slices
% 1 - Sagittal, 2 - Coronal, 3 - Axial
Plane2d = 3;     

% If true, crops the 2d image so that the FOV becomes CropSize2d smalled in
% each direction
Crop2d     = false;
CropSize2d = 20;

% If true, pull out a slab from the 3D image of size 2*SlabSize3d
ExtractSlab3d = false;
SlabSize3d    = 5;

% % Define thick-slices
% DownSampling = 1/6;
% TS           = {[DownSampling 1 1], ... 
%                 [1 DownSampling 1], ...
%                 [2*DownSampling 2*DownSampling 1;]};
% Gap          = 0;

%--------------------------------------------------------------------------
% Create output directory
%--------------------------------------------------------------------------

DirSim3D = fullfile(DirSim,'3D');
DirSim2D = fullfile(DirSim,'2D');
if  (exist(DirSim3D,'dir') == 7),  rmdir(DirSim3D,'s'); end; mkdir(DirSim3D);
if  (exist(DirSim2D,'dir') == 7),  rmdir(DirSim2D,'s'); end; mkdir(DirSim2D);

if ~isempty(DirPrint) 
    if (exist(DirPrint,'dir') == 7)  
        rmdir(DirPrint,'s'); 
    end
    mkdir(DirPrint);
end

%--------------------------------------------------------------------------
% Get reference BrainWeb images
%--------------------------------------------------------------------------

if isempty(DirRef)
    Nii_ref = nifti(spm_select(Inf,'nifti','Select image')); 
else
    Nii_ref = nifti(spm_select('FPList',DirRef,'^.*\.nii$'));
end
C = numel(Nii_ref); % Number of channels


if isempty(Offset),       Offset       = repmat({[0; 0; 0]},1,C); end
if isempty(Rotation),     Rotation     = repmat({[0; 0; 0]},1,C); end
if isempty(DownSampling), DownSampling = repmat({[0 0 0]},1,C); end
if isempty(mrg),          mrg          = repmat({[0 0 0]},1,C); end

%--------------------------------------------------------------------------
% Simulate
%--------------------------------------------------------------------------

 % For storing filenames
fnames_in    = cell(1,C);
fnames_out3d = cell(1,C);
fnames_out2d = cell(1,C);
bb0          = repmat({[1 1 1; 1 1 1]},1,C);
for c=1:C % Loop over channels

    fname       = Nii_ref(c).dat.fname;    
    [~,nam,ext] = fileparts(fname);  
    
    fnames_in{c} = fname;
    
    mat = Nii_ref(c).mat;
    img = Nii_ref(c).dat(:,:,:);   
    dm  = size(img);
    
    if BiasFieldScl
        % Multiply with bias field        
        vx  = sqrt(sum(mat(1:3,1:3).^2));
        dm  = size(img);
        bf  = sample_bf(BiasFieldScl,dm,vx);
        img = bf.*img;
        clear bf
    end
    
    if any(DownSampling{c} > 0) && any(DownSampling{c} ~= 1)
        % Down-sample image w NN interpolation
        [img,mat] = resample_img(Nii_ref(c),img,DownSampling{c},DownSamplingDeg);      
        dm        = size(img);
    end
        
%     if ThickSliced
%         % Build dat object       
%         D    = diag([TS{c} 1]);
%         mat1 = mat/D;
%         dm1  = floor(D(1:3,1:3)*dm')';
%         vx   = sqrt(sum(mat1(1:3,1:3).^2));   
% 
%         proj.mat     = mat1;
%         proj.dat.dim = dm1;
% 
%         dat = init_dat('superres',{proj},mat,dm,[],Gap);
%         img = A(single(img),dat);
%         img = img{1};
%         mat = mat1;
%     end            

    if NoisePrct > 0
        % Add noise 
        % MR images are usually reconstructed as the magnitude of
        % an image that was originally complex.  Gaussian noise on complex
        % data gives a Rice distribution in the magnitude.
        
        msk = isfinite(img) & img > 0;
        mn  = mean(img(msk));
        
%         img = abs(img + (NoisePrct*mn)*randn(size(img)));        
%         img = abs(img + (NoisePrct*mn)*randn(size(img)) + sqrt(-1)*(NoisePrct*mn)*randn(size(img)));
        v1  = 0;
        s1  = NoisePrct*mn;
        img = img + ricernd(v1*ones(size(img)), s1);
    end      
    
    %----------------------------------------------------------------------
    % Write to NIfTI
    %----------------------------------------------------------------------
        
    % 3D
    nfname = fullfile(DirSim3D,[nam ext]); 
    create_nii(nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'Simulated (3D)');

    fnames_out3d{c} = nfname;
end % End loop over channels

for c=1:C % Loop over channels
        
    % 2D
    fname2d         = extract_slice(fnames_out3d{c},Plane2d);
    [~,nam,ext]     = fileparts(fname2d);  
    movefile(fname2d,DirSim2D);
    fname2d         = fullfile(DirSim2D,[nam ext]);    
    fnames_out2d{c} = fname2d;
    
    if any(mrg{c} ~= 0)
        % 3D
        ofname = fnames_out3d{c};
        Nii    = nifti(ofname);
        mat    = Nii.mat;      
        vx     = sqrt(sum(mat(1:3,1:3).^2));
        
        mrg{c} = mrg{c}./vx;
        
        dm              = Nii.dat.dim;      
        bb              = [1 + mrg{c}(1) dm(1) - mrg{c}(1); 1 + mrg{c}(2) dm(2) - mrg{c}(2); 1 + mrg{c}(3) dm(3) - mrg{c}(3);]';
        [nfname,bb]     = subvol(spm_vol(ofname),bb);        
        fnames_out3d{c} = nfname;    
        delete(ofname);        
        
        bb0{c} = bb;

%         Nii    = nifti(nfname);
%         mat    = Nii.mat; 
        
        % 2D
        ofname2d        = fnames_out2d{c};
        Nii             = nifti(ofname2d);
        dm              = Nii.dat.dim;        
        bb              = [1 + mrg{c}(1) dm(1) - mrg{c}(1); 1 + mrg{c}(2) dm(2) - mrg{c}(2); -Inf Inf;]';
        fname2d         = subvol(spm_vol(ofname2d),bb);
        fnames_out2d{c} = fname2d;
        delete(ofname2d);        
    end

end % End loop

R3d = repmat(eye(4),[1 1 C]);
R2d = R3d;
for c=1:C
    R3d(:,:,c) = rigidly_realign(fnames_out3d{c},Offset{c},Rotation{c});
    R2d(:,:,c) = rigidly_realign(fnames_out2d{c},Offset{c},Rotation{c});
end

if Reslice
    fnames_out3d = reslice_imgs(fnames_out3d,Padding);
end

oNii   = nifti;
oNii2d = nifti;
for c=1:C
    oNii(c)   = nifti(fnames_out3d{c});
    oNii2d(c) = nifti(fnames_out2d{c});
end

if ShowSimulated   
    % Show simulated results

    fnames = {};
    for c=1:C
%         fnames{end + 1} = fnames_in{c};
        fnames{end + 1} = fnames_out3d{c};
    end

    spm_check_registration(char(fnames))
end
%==========================================================================

%==========================================================================
function fname = extract_slice(fname,ax)
% Extract a 2D slice from a 3D volume. The variable ax (1,2,3) decides from
% what dimension to pull out the slice.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   

if nargin < 2, ax = 3; end

% Create bounding box
V  = spm_vol(fname);
dm = V.dim;
if ax     == 1
    d1 = floor(dm(1)/2) + 1;
    bb = [d1 d1;-inf inf;-inf inf];   
elseif ax == 2
    d1 = floor(dm(2)/2) + 1;
    bb = [-inf inf;d1 d1;-inf inf];
elseif ax == 3 
    d1 = floor(dm(3)/2) + 1;
    bb = [-inf inf;-inf inf;d1 d1];
end                

% Crop according to bounding-box
fname = subvol(V,bb','2d_');

% Make sure 1D plane is in z dimension
Nii  = nifti(fname);
mat  = Nii.mat;
img  = Nii.dat(:,:,:);

if ax == 1 || ax == 2
    % Permute image data and apply permutation matrix to orientation matrix
    if ax == 1
        img = permute(img,[2 3 1]);            
        P   = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    else
        img = permute(img,[1 3 2]);        
        P   = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    end   

    mat     = P*mat*P';    
end

% Set zero translation in z-direction
mat(3,4) = 0;

% Overwrite image data
VO             = spm_vol(fname);
dm             = [size(img) 1];
VO.dim(1:3)    = dm(1:3);        
VO.mat         = mat;
VO             = spm_create_vol(VO);        
Nii            = nifti(VO.fname);    
Nii.dat(:,:,:) = img; 
%==========================================================================

%==========================================================================
function R = rigidly_realign(fname,offset,rotation)
% Rigidly (translation and rotation) modify an orientation matrix
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Nii  = nifti(fname);
dm   = Nii.dat.dim;
is2d = numel(dm) == 2;

if is2d
    q = [offset(1) offset(2) 0 0 0 rotation(3)];
else
    q = [offset(1) offset(2) offset(3) rotation(1) rotation(2) rotation(3)];
end
R  = spm_matrix(q(:)');
MM = spm_get_space(fname);
spm_get_space(fname, R*MM);
%==========================================================================

%==========================================================================
function bf = sample_bf(scl,lat,vs,fwhm,prm,Verbose)
% Sample a bias-field.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 1, scl      = 1; end
if nargin < 2, lat      = [256 256 150]; end
if nargin < 3, vs       = [1 1 1]; end
if nargin < 4, fwhm     = 100; end
if nargin < 5, prm      = [1e-3 0 1e-6]; end
if nargin < 6, Verbose  = false; end

if scl == 0
    bf = ones(lat);
    return
end

nbcmp = fwhm2nbcomp(lat, vs, fwhm);

L     = regulariser(prm, lat, nbcmp, vs, 'n');
S     = inv(L);
[U,S] = svd(S);
S     = sqrt(S);
c     = reshape(U*S*randn(size(L,1),1), nbcmp);

[B1,B2,B3] = dcbasis(lat, nbcmp);

bf = zeros(lat);
for z=1:lat(3)
    bf(:,:,z) = exp(scl*transf(B1,B2,B3(z,:),c));
end

if Verbose || nargin == 0 || nargin == 1
    figure(111)
    imagesc3d(bf); colorbar
end
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%==========================================================================

%==========================================================================
function nbcmp = fwhm2nbcomp(lattice, vs, fwhm)
% FORMAT nbcmp = spm_bias_lib('fwhm2nbcomp', lattice, voxel_size, fwhm)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
% fwhm         - Full-width half-max of the highest frequency basis (mm)
%
% The number of components is chosen so that the full-width half-max of 
% the highest frequency basis function is smaller than fwhm. The bias 
% field cannot model effects whose spatial frequency is higher than this 
% value.
%
% If only one value is provided for voxel_size or fwhm, the same value is
% used along all dimensions.

% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end
fwhm = reshape(fwhm, 1, []);
if numel(fwhm) < ndim
    fwhm = padarray(fwhm, [0 ndim-numel(fwhm)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute number of components per direction
nbcmp = ceil(2 * vs .* lattice ./ fwhm);
nbcmp = max(nbcmp, 1);
%==========================================================================

%==========================================================================
function varargout = dcbasis(lattice, nb_component)
% FORMAT [Bx,By,Bz,...] = spm_bias_lib('dcbasis', lattice, nb_component)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
%
% Bx - Smooth basis along the x dimension [dx*nx] 
% By - Smooth basis along the y dimension [dy*ny]
% ...
%
% There are as many basis objects as elements in `lattice`

ndim = numel(lattice);

% -------------------------------------------------------------------------
% Preprocess input arguments
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute each basis
varargout = cell(1,min(ndim, nargout));
for d=1:min(ndim, nargout)
    varargout{d} = spm_dctmtx(lattice(d),nb_component(d));
end
%==========================================================================

%==========================================================================
function L = regulariser(mode, lattice, nb_component, vs, bnd)
% FORMAT L = regulariser(param, lattice, nb_component, voxel_size)
% FORMAT L = regulariser(mode,  lattice, nb_component, voxel_size)
%
% param        - Parameters for absolute, membrane and bending energies
% mode         - Name of a single energy ('absolute'/'membrane'/'bending')
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
%
% L            - Precision matrix [(nx*ny*...)^2]
%
% If numerical parameters are provided, a weighted combination of the  
% three types of regularisation is returned.
% If an energy name is provided, the matrix that allows to compute it is
% returned (without weighting: the regularisation parameter should be 
% multiplied with this matrix)
%
% If only one value is provided for nb_component or voxel_size, the
% same value is used along all dimensions.

if nargin < 5
    bnd = 'neumann';
end

% -------------------------------------------------------------------------
% Special case: mixture of regularisers
if ~ischar(mode)
    param = mode;
    L = 0;
    for i=1:numel(param)
        if param(i) ~= 0
            switch i
                case 1
                    mode = 'absolute';
                case 2
                    mode = 'membrane';
                case 3
                    mode = 'bending';
                case 4
                    mode = 'linearelastic1';
                case 5
                    mode = 'linearelastic2';
            end
            L1 = param(i) * regulariser(mode, lattice, nb_component, vs, bnd);
            if numel(L) == 1 || size(L,1) == size(L1,1)
                L = L + L1;
            else
                L0   = L;
                nprm = size(L,1);
                ndim = size(L1,1)/nprm;
                L = zeros(ndim*nprm);
                for d=1:ndim
                    L(((d-1)*nprm+1):d*nprm,((d-1)*nprm+1):d*nprm) = L0;
                end
                clear L0
                L = L + L1;
            end
        end
    end
    return
end


% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end
if nargin < 4
    vs = 1;
end
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Mode-specific options
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        maxdiff = 0;
    case {'membrane' 'mem' 'm' ...
          'linear-elastic1' 'linearelastic1' 'le1' ...
          'linear-elastic2' 'linearelastic2' 'le2'}
        maxdiff = 1;
    case {'bending' 'ben' 'b'}
        maxdiff = 2;
    otherwise
        error('Unknown mode %s, should be ''absolute'', ''membrane'' or ''bending''.', mode);
end

% -------------------------------------------------------------------------
% Compute each basis + square it
switch lower(bnd)
    case {0, 'circulant', 'circ', 'c'}
        mtxfun = @spm_dftmtx;
    case {1, 'neumann', 'neu', 'n'}
        mtxfun = @spm_dctmtx;
    case {2, 'dirichlet', 'dir', 'd'}
        mtxfun = @spm_dstmtx;
    otherwise
        error('Unknown boundary condition');
end

basis = cell(ndim, maxdiff + 1);
nbprm = 1;
for d=1:ndim
    for diff=0:maxdiff
        switch diff
            case 0
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d));
                nbprm = nbprm * size(basis{d,diff+1}, 2);
            case 1
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff') / vs(d);
            case 2
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff2') / vs(d)^2;
        end
        if any(strcmpi(mode, {'absolute' 'abs' 'a' 'membrane' 'mem' 'm' 'bending' 'ben' 'b'}))
            basis{d,diff+1} = basis{d,diff+1}.' * basis{d,diff+1};
        end
    end
end

% -------------------------------------------------------------------------
% Compute precision matrix
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        L = 1;
        for d=1:ndim
            L = spm_krutil(basis{d,1}, L);
        end
    case {'membrane' 'mem' 'm'}
        L = 0;
        for dd=1:ndim               % Which dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd
                    L1 = spm_krutil(basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
        end
    case {'bending' 'ben' 'b'}
        L = 0;
        for dd1=1:ndim              % First dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd1
                    L1 = spm_krutil(basis{d,3}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
            for dd2=dd1+1:ndim      % Second dimension to differentiate
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd1 || d == dd2
                        L1 = spm_krutil(basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}, L1);
                    end
                end
                L = L + 2 * L1;
            end
        end
    case {'linear-elastic1' 'linearelastic1' 'le1'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            for dd=1:ndim          % First dimension to differentiate
                if dd == h1
                    coeff = 1;
                else
                    coeff = 0.5;
                end
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd
                        L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h1) = L(:,h1,:,h1) + coeff * reshape(L1, [nbprm 1 nbprm]);
            end
            
            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end
            
        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
    case {'linear-elastic2' 'linearelastic2' 'le2'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            L1 = 1;
            for d=1:ndim        % Kronecker loop
                if d == h1
                    L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                end
            end
            L(:,h1,:,h1) = L(:,h1,:,h1) + 0.5 * reshape(L1, [nbprm 1 nbprm]);
            
            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end
            
        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
end
%==========================================================================

%==========================================================================
function [fname,bb] = subvol(V,bb,prefix)
% Extract a subvolume
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   

if nargin < 3, prefix = 'sv_'; end

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
bb(2,:) = min(bb(2,:),V.dim(1:3));

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
fname         = fullfile(pth,[prefix nam ext]);
VO.fname      = fname;
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3)
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),0);
    VO  = spm_write_plane(VO,img,z);
end
%==========================================================================

%==========================================================================
function r = ricernd(v, s)
%RICERND Random samples from the Rice/Rician probability distribution.
%   r = ricernd(v, s) returns random sample(s) from the Rice (aka Rician) 
%   distribution with parameters v and s.
%   (either v or s may be arrays, if both are, they must match in size)
%
%   R ~ Rice(v, s) if R = sqrt(X^2 + Y^2), where X ~ N(v*cos(a), s^2) and
%   Y ~ N(v*sin(a), s^2) are independent normal distributions (any real a).
%   Note that v and s are *not* the mean and standard deviation of R!
%
%   The size of Y is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   Note, to add Rician noise to data, with given s and data-dependent v:
%     new = ricernd(old, s);
%
%   Reference: http://en.wikipedia.org/wiki/Rice_distribution (!)
%
%   Example:
%
%     % Compare histogram of random samples with theoretical PDF:
%     v = 4; s = 3; N = 1000;
%     r = ricernd(v*ones(1, N), s);
%     c = linspace(0, ceil(max(r)), 20);
%     w = c(2); % histogram bin-width
%     h = histc(r, c); bar(c, h, 'histc'); hold on
%     xl = xlim; x = linspace(xl(1), xl(2), 100);
%     plot(x, N*w*ricepdf(x, v, s), 'r');
%     
%   See also RICEPDF, RICESTAT, RICEFIT

%   Missing (?) 'See also's RICECDF, RICEINV, RICELIKE

%   Inspired by normpdf from the MATLAB statistics toolbox
%   Copyright 2008 Ged Ridgway (Ged at cantab dot net)

if isscalar(v)
    dim = size(s);
elseif isscalar(s)
    dim = size(v);
elseif all(isequal(size(v), size(s)))
    % (both non-scalar, matching)
    dim = size(v); % == size(s)
else
    error('ricernd:InputSizeMismatch','Sizes of s and v inconsistent.')
end

x = s .* randn(dim) + v;
y = s .* randn(dim);
r = sqrt(x.^2 + y.^2);
%==========================================================================

%==========================================================================
function [img,mat,dm] = resample_img(Nii,img,samp,deg,bc)
% Resample an image using deg interpolation, with bc boundary conditions.
% If samp < 1, does down-sampling; if samp > 1, does up-sampling.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 4, deg = 0; end
if nargin < 5, bc  = 0; end

if numel(samp) == 1, samp = samp*ones([1 3]); end
if numel(deg)  == 1, deg  = deg*ones([1 3]);  end
if numel(bc)   == 1, bc   = bc*ones([1 3]);   end

samp(samp == 0) = 1;

if numel(Nii.dat.dim) ~= 3
    error('numel(Nii.dat.dim) ~= 3')
end

% Input image properties
mat0 = Nii.mat;
dm0  = size(img);

% Output image properties
D    = diag([samp 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0')';

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

% Resample
img                         = spm_bsplins(img,x1,y1,z1,[deg bc]);    
img(~isfinite(img) | img<0) = 0;
%==========================================================================

%==========================================================================
function Nii = create_nii(pth,dat,mat,dtype,descrip,offset,scl_slope,scl_inter)
% Create a NIfTI file
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

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

Nii.dat(:) = dat(:);
%==========================================================================
