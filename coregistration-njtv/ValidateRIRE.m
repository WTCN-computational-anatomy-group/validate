function e = ValidateRIRE(DirData,DirRes)
if nargin < 1, DirData = 'test-data/RIRE/training_001'; end
if nargin < 2, DirRes  = 'test-data/RIRE/training_001/results'; end
    
DO_WRITE  = false;
DirTransf = fullfile(DirData,'transformations');

% Get images
[Nii,info] = LoadRIRE(DirData);

% Load rigid matrices
PthR = fullfile(DirRes,'R.mat');
var  = load(PthR);
R    = var.R0;

% Select from and to scans
from = struct('ix',0,'nam','');
to   = from;
fn   = fieldnames(info);
cf   = 1;
ct   = 1;
for i=1:numel(fn)
    if strcmp(fn{i},'ct') || strcmp(fn{i},'pet')
        from(cf) = info.(fn{i});
        cf       = cf + 1;
    else
        to(ct) = info.(fn{i});
        ct     = ct + 1;
    end
end

% map_aec        = containers.Map;
% map_aec('ct')  = [];
% map_aec('pet') = [];
for i1=1:numel(from) % loop over from images
    
    % Get from indices and names        
    ix_from  = from(i1).ix;
    nam_from = from(i1).nam;
        
%     aec = map_aec(nam_from);
    
    for i2=1:numel(to) % loop over to images
        
        % Get to indices and names        
        ix_to  = to(i2).ix;        
        nam_to = to(i2).nam;

        % Get RIRE transformation files        
        PthTransf = GetPthTransf(DirTransf,nam_from,nam_to);

        % Compute error (when no rigid registration has been done)
        [~,~,ecid] = ComputeError(Nii(ix_from),eye(4),eye(4),PthTransf);
        
        % Get from and to rigid transformations
        R_from = R(:,:,ix_from);
        R_to   = R(:,:,ix_to);

        % Compute error
        [c_from,c_to,ec] = ComputeError(Nii(ix_from),R_from,R_to,PthTransf);
          
        e         = struct;
        e.from_to = [nam_from '-' nam_to];        
        e.median  = [median(ecid) median(ec)];
        e.max     = [max(ecid) max(ec)];
%         e.mean   = mean(ec);        
%         e.c       = ec;
%         e.N      = numel(ec);

        % show error
        disp(e);

        if DO_WRITE
            % write output
            WriteRes(DirRes,c_from,c_to,nam_from,nam_to);
        end
    end    
%     map_aec(nam_from) = aec;
end
%==========================================================================

%==========================================================================
function [c_from_o,c_to,ec] = ComputeError(Nii_from,R_from,R_to,PthTransf)

GroundTruthGiven = false;
if isfile(PthTransf), GroundTruthGiven = true; end

if GroundTruthGiven
    % Get ground-truth
    [c_from0,c_to0] = GetGroundTruth(PthTransf);
end

% Parameters
M_from = Nii_from.mat;
d_from = Nii_from.dat.dim;
% v_from = sqrt(sum(MatFrom(1:3,1:3).^2));

M_rire = GetRIRESpace;

% CornersFromEst
c_from = [1    1    1    1
          1    1    d_from(3) 1
          1    d_from(2) 1    1
          1    d_from(2) d_from(3) 1
          d_from(1) 1    1    1
          d_from(1) 1    d_from(3) 1
          d_from(1) d_from(2) 1    1
          d_from(1) d_from(2) d_from(3) 1]'; 

ix     = [1 5 3 7 2 6 4 8];
c_from = c_from(:,ix);

c_from_o = M_from*c_from;
c_from_o = M_rire*c_from_o;
c_from_o = c_from_o(1:3,:);

if GroundTruthGiven && 0, disp(c_from_o - c_from0); end % correct?

c_from = M_from*c_from;

R = R_to\R_from;
% R0 = GetMatTrue(CornersFromTrue,CornersToTrue);

c_to = R\c_from;
c_to = M_rire*c_to;
c_to = c_to(1:3,:);

if GroundTruthGiven
    % Corner errors
    ec = abs(c_to - c_to0);
    ec = ec(:);
else
    ec = 0;
end

c_from_o = round(c_from_o,4);
c_to     = round(c_to,4);
%==========================================================================

%==========================================================================
function WriteRes(DirRes,c_from,c_to,nam_from,nam_to)

PthRes = fullfile(DirRes,[nam_from '_' nam_to '.standard']);
if isfile(PthRes), delete(PthRes); end

% Get template
Text = GetTemplate;

if ~contains(nam_from,'ct') && ~contains(nam_from,'pet')
    nam_from = ['mr_' nam_from];
end
if ~contains(nam_to,'ct') && ~contains(nam_to,'pet')
    nam_to = ['mr_' nam_to];
end

% Date
ix       = 9;
t        = string(Text{ix}(1:6));
s        = date;
t        = strcat(t,s);
Text{ix} = char(t);

% Patient number
num_pat = strsplit(DirRes,filesep);
num_pat = num_pat{end - 1};
num_pat = strsplit(num_pat,'_');
num_pat = num_pat{end};

ix       = 10;
t        = string(Text{ix}(1:16));
s        = num_pat;
t        = strcat(t,s);
Text{ix} = char(t);

% From
ix       = 11;
t        = string(Text{ix}(1:6));
s        = nam_from;
t        = strcat(t,s);
Text{ix} = char(t);

% To
ix       = 12;
t        = string(Text{ix}(1:4));
s        = nam_to;
t        = strcat(t,s);
Text{ix} = char(t);

% Coordinates
for c=1:8
    ix       = 15 + c;
    t        = string(Text{ix}(1:10));
    s        = [c_from(1:3,c)' c_to(1:3,c)'];
    s        = mat2str(s);
    s        = s(2:end - 1);
    t        = strcat(t,s);
    Text{ix} = char(t);
end

% Write
fid = fopen(PthRes, 'w');
for i = 1:numel(Text)
    if Text{i+1} == -1
        fprintf(fid,'%s', Text{i});
        break
    else
        fprintf(fid,'%s\n', Text{i});
    end
end
%==========================================================================

%==========================================================================
function Text = GetTemplate
PthTemplate = 'test-data/RIRE/template.standard';

% Read txt into cell A
fid = fopen(PthTemplate, 'r');
i = 1;
tline = fgetl(fid);
Text{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    Text{i} = tline;
end
fclose(fid);
%==========================================================================

%==========================================================================
function M_rire = GetRIRESpace
% X = CornersFromTrue;
% Y = CornersFromEst;
% X = [X; ones(1,8)];
% Y = [Y; ones(1,8)];
% T = Y/X;
M_rire = [-1  0 0 0; ...
           0 -1 0 0; ...
           0  0 1 0; ...
           0  0 0 1];  
%==========================================================================

%==========================================================================
function [c_from0,c_to0] = GetGroundTruth(PathTransf)
fid  = fopen(PathTransf,'rt');
Text = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

Text            = Text{1};
c_from0 = zeros(3,8);
c_to0   = zeros(3,8);
for i=1:8
    vals = strsplit(Text{15 + i});
    
    c_from0(1,i) = str2double(vals{2});
    c_from0(2,i) = str2double(vals{3});
    c_from0(3,i) = str2double(vals{4});
    
    c_to0(1,i) = str2double(vals{5});
    c_to0(2,i) = str2double(vals{6});
    c_to0(3,i) = str2double(vals{7});
end
%==========================================================================

%==========================================================================
function PthTransf = GetPthTransf(DirTransf,nam_from,nam_to)
PthTransf = fullfile(DirTransf,[nam_from '_' nam_to '.standard']);
%==========================================================================

% %==========================================================================
% function MatTrue = GetMatTrue(CornersFromTrue,CornersFromEst)
% X = CornersFromTrue;
% Y = CornersFromEst;
% X = [X; ones(1,8)];
% Y = [Y; ones(1,8)];
% MatTrue = Y/X;
% %==========================================================================