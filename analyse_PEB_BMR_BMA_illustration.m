% load and examine example of a fully connected, not inverted code
load('DCM_example_not_inverted.mat')

%% copy it across participants
% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'/Users/user/specify_group_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
%% invert the fully connected model

for s = 1:size(subcode,1)
clc
s
m1 = cellstr([root,filesep,subcode(s,:), filesep, 'results', filesep, 'DCM_full_m0001.mat']);


clear matlabbatch

matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = m1;

spm_jobman('run',matlabbatch);

end

disp('inverted model for all subjects')
%% grab GCM_HC and GCM_NF1, make GLM, analyse A matrix


load('/Users/user/data/GCM_full.mat', 'GCM')

M   = struct();
M.Q = 'all'; 
N = 16 + 28; % always conrols first
M.X(:,1) = ones(N, 1);
M.X(:,2) = ones(N,1);
M.X(1:16,2) = -1; % note: controls have -1 value, patients have 1

load('demographics.mat')

age = [age_HC; age_NF1];
sex = [sex_HC; sex_NF1];

M.X(:,3) = age; % age at first visit
M.X(:,4) = sex; % 1 = male, 0 = female

% Choose field
field = {'A'};

% Estimate model
PEB     = spm_dcm_peb(GCM,M,field);

save('PEB_A.mat','PEB');



%%  Switch off model connections for PEB
load('/Users/user/GCM_HC_NF1.mat')
for s = 1:N % switch off model connections and place these templates next to each other (sub:DCM)
% Get an existing model. We'll use the first subject's DCM as a template
DCM_full = load(GCM{s});

% IMPORTANT: Clear out the old priors, or changes to DCM.a,b,c will be ignored
if isfield(DCM_full,'M')
    DCM_full = rmfield(DCM_full ,'M');
end

% Specify candidate models that differ in particular A-matrix connections, e.g.
m1 = DCM_full;
m1.b(1,1) = 0; % Switching off the connection from region 1 to region 1
m2 = DCM_full;
m2.b(1,2) = 0; % Switching off the connection from region 2 to region 1
m3 = DCM_full;
m3.b(1,3) = 0; 
m4 = DCM_full;
m4.b(1,4) = 0; 
m5 = DCM_full;
m5.b(1,5) = 0; 
m6 = DCM_full;
m6.b(1,6) = 0; 

m7 = DCM_full;
m7.b(2,1) = 0; 
m8 = DCM_full;
m8.b(2,2) = 0; 
m9 = DCM_full;
m9.b(2,3) = 0; 
m10 = DCM_full;
m10.b(2,4) = 0; 
m11 = DCM_full;
m11.b(2,5) = 0; 
m12 = DCM_full;
m12.b(2,6) = 0; 

m13 = DCM_full;
m13.b(3,1) = 0; 
m14 = DCM_full;
m14.b(3,2) = 0; 
m15 = DCM_full;
m15.b(3,3) = 0; 
m16 = DCM_full;
m16.b(3,4) = 0; 
m17 = DCM_full;
m17.b(3,5) = 0; 
m18 = DCM_full;
m18.b(3,6) = 0; 

m19 = DCM_full;
m19.b(4,1) = 0; 
m20 = DCM_full;
m20.b(4,2) = 0; 
m21 = DCM_full;
m21.b(4,3) = 0; 
m22 = DCM_full;
m22.b(4,4) = 0; 
m23 = DCM_full;
m23.b(4,5) = 0; 
m24 = DCM_full;
m24.b(4,6) = 0; 

m25 = DCM_full;
m25.b(5,1) = 0; 
m26 = DCM_full;
m26.b(5,2) = 0; 
m27 = DCM_full;
m27.b(5,3) = 0; 
m28 = DCM_full;
m28.b(5,4) = 0; 
m29 = DCM_full;
m29.b(5,5) = 0; 
m30 = DCM_full;
m30.b(5,6) = 0; 

m1 = DCM_full;
m1.b(6,1) = 0; 
m2 = DCM_full;
m2.b(6,2) = 0; 
m3 = DCM_full;
m3.b(6,3) = 0; 
m4 = DCM_full;
m4.b(6,4) = 0; 
m5 = DCM_full;
m5.b(6,5) = 0; 
m6 = DCM_full;
m6.b(6,6) = 0; 

GCM_out(s,:) =  {DCM_full, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18 ...
     m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29,m30};

end

clear GCM
GCM = GCM_out;
save('GCM_templates.mat', 'GCM')
%% Search over nested PEB models
BMA = spm_dcm_peb_bmc(PEB);
spm_dcm_peb_review(BMA,GCM) % results of BMA
save('BMA_GCM_A.mat', 'GCM', 'BMA')
%%
% safe_BMA = BMA; % in case you need to explore the results and do not want
% to load BMA file 

BMA.Ep(BMA.Pp<.95) = 0;
pos = reshape(BMA.Ep(37:72), 6, 6);
pos = full(pos);
pos = round(pos, 2);
pos(pos<0) = 0;
A = ones(6);
A(pos == 0) = 0;
G = digraph(A, names);

% reordering everything, for intutive graph
order = {'right Precuneus' 'right dlPFC' 'left dlPFC'  'left SMA' 'left IPG' 'right IPG' };
G = reordernodes(G,order);


pos = digraph(pos, names);
pos = reordernodes(pos,order);

LWidths = 5*pos.Edges.Weight/max(pos.Edges.Weight); % if you want to 

EL = pos.Edges.Weight ; % this will be replaced with model number


figure(1)
p = plot(G,'Layout','circle', ...
    'NodeFontSize', 15, ...
    'EdgeLabel',EL, 'EdgeLabelColor', 'g', 'EdgeFontSize', 12)

p.EdgeColor = 'g';
p.ArrowSize=10;
p.NodeColor = 'k';
p.ArrowPosition = 0.99;

%% illustrate only negative connections
names = {'left IPG' 'left SMA' 'left dlPFC' 'right IPG', 'right Precuneus', 'right dlPFC'};

BMA.Ep(BMA.Pp<.95) = 0;
pos = reshape(BMA.Ep(37:72), 6, 6);
pos = full(pos);
pos = round(pos, 2);
pos = pos
pos(pos>0) = 0;

A = ones(6);
A(pos == 0) = 0;
G = digraph(A, names);

% reordering everything, for intutive graph
order = {'right Precuneus' 'right dlPFC' 'left dlPFC'  'left SMA' 'left IPG' 'right IPG' }
G = reordernodes(G,order);


pos = digraph(pos, names);
pos = reordernodes(pos,order);

LWidths = 5*pos.Edges.Weight/max(pos.Edges.Weight);

EL = pos.Edges.Weight ; % this will be replaced with model number



figure(2)
p = plot(G,'Layout','circle', ...
    'NodeFontSize', 15, ...
    'EdgeLabel',EL, 'EdgeLabelColor', 'r', 'EdgeFontSize', 12)
p.EdgeColor = 'r';
p.ArrowSize=10;
p.NodeColor = 'k';
p.ArrowPosition = 0.99;

%% repeat the lines 51:235 for B matix
load('/Users/user/data/GCM_full.mat', 'GCM')

M   = struct();
M.Q = 'all'; 
N = 16 + 28; % always conrols first
M.X(:,1) = ones(N, 1);
M.X(:,2) = ones(N,1);
M.X(1:16,2) = -1; % note: controls have -1 value, patients have 1

load('demographics.mat')

age = [age_HC; age_NF1];
sex = [sex_HC; sex_NF1];

M.X(:,3) = age; % age at first visit
M.X(:,4) = sex; % 1 = male, 0 = female

% Choose field
field = {'B'};

% Estimate model
PEB     = spm_dcm_peb(GCM,M,field);

save('PEB_B.mat','PEB');



%%  Switch off model connections for PEB
load('/Users/user/GCM_HC_NF1.mat')
for s = 1:N % switch off model connections and place these templates next to each other (sub:DCM)
% Get an existing model. We'll use the first subject's DCM as a template
DCM_full = load(GCM{s});

% IMPORTANT: Clear out the old priors, or changes to DCM.a,b,c will be ignored
if isfield(DCM_full,'M')
    DCM_full = rmfield(DCM_full ,'M');
end

% Specify candidate models that differ in particular A-matrix connections, e.g.
m1 = DCM_full;
m1.b(1,1) = 0; % Switching off the connection from region 1 to region 1
m2 = DCM_full;
m2.b(1,2) = 0; % Switching off the connection from region 2 to region 1
m3 = DCM_full;
m3.b(1,3) = 0; 
m4 = DCM_full;
m4.b(1,4) = 0; 
m5 = DCM_full;
m5.b(1,5) = 0; 
m6 = DCM_full;
m6.b(1,6) = 0; 

m7 = DCM_full;
m7.b(2,1) = 0; 
m8 = DCM_full;
m8.b(2,2) = 0; 
m9 = DCM_full;
m9.b(2,3) = 0; 
m10 = DCM_full;
m10.b(2,4) = 0; 
m11 = DCM_full;
m11.b(2,5) = 0; 
m12 = DCM_full;
m12.b(2,6) = 0; 

m13 = DCM_full;
m13.b(3,1) = 0; 
m14 = DCM_full;
m14.b(3,2) = 0; 
m15 = DCM_full;
m15.b(3,3) = 0; 
m16 = DCM_full;
m16.b(3,4) = 0; 
m17 = DCM_full;
m17.b(3,5) = 0; 
m18 = DCM_full;
m18.b(3,6) = 0; 

m19 = DCM_full;
m19.b(4,1) = 0; 
m20 = DCM_full;
m20.b(4,2) = 0; 
m21 = DCM_full;
m21.b(4,3) = 0; 
m22 = DCM_full;
m22.b(4,4) = 0; 
m23 = DCM_full;
m23.b(4,5) = 0; 
m24 = DCM_full;
m24.b(4,6) = 0; 

m25 = DCM_full;
m25.b(5,1) = 0; 
m26 = DCM_full;
m26.b(5,2) = 0; 
m27 = DCM_full;
m27.b(5,3) = 0; 
m28 = DCM_full;
m28.b(5,4) = 0; 
m29 = DCM_full;
m29.b(5,5) = 0; 
m30 = DCM_full;
m30.b(5,6) = 0; 

m1 = DCM_full;
m1.b(6,1) = 0; 
m2 = DCM_full;
m2.b(6,2) = 0; 
m3 = DCM_full;
m3.b(6,3) = 0; 
m4 = DCM_full;
m4.b(6,4) = 0; 
m5 = DCM_full;
m5.b(6,5) = 0; 
m6 = DCM_full;
m6.b(6,6) = 0; 

GCM_out(s,:) =  {DCM_full, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18 ...
     m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29,m30};

end

clear GCM
GCM = GCM_out;
save('GCM_templates.mat', 'GCM')
%% Search over nested PEB models
BMA = spm_dcm_peb_bmc(PEB);
spm_dcm_peb_review(BMA,GCM) % results of BMA
save('BMA_GCM_B.mat', 'GCM', 'BMA')
%%
% safe_BMA = BMA; % in case you need to explore the results and do not want
% to load BMA file 

BMA.Ep(BMA.Pp<.95) = 0; % threshold connections by their posterior probability
pos = reshape(BMA.Ep(37:72), 6, 6); % only take Covariate 2, which is the effect of diagnosis
pos = full(pos); % sparse matrix to full matrix
pos = round(pos, 2); % rounding it for clarity
pos(pos<0) = 0; % only illustrate the positive connections
A = ones(6); % A-matrix
A(pos == 0) = 0; % binarise A matrix
G = digraph(A, names); % for makes a binary digraph

% reordering everything, for a visually more intutive graph
order = {'right Precuneus' 'right dlPFC' 'left dlPFC'  'left SMA' 'left IPG' 'right IPG' };
G = reordernodes(G,order);


pos = digraph(pos, names); % reordering also the labels
pos = reordernodes(pos,order);

LWidths = 5*pos.Edges.Weight/max(pos.Edges.Weight); % if you want to value line width according to its weight

EL = pos.Edges.Weight ;


figure(1)
p = plot(G,'Layout','circle', ...
    'NodeFontSize', 15, nodeColor, 'k', ...
    'EdgeLabel',EL, 'EdgeLabelColor', 'g', 'EdgeFontSize', 12, ...
    'ArrowSize', 10, 'ArrowPosition', 0.99)

%% illustrate only negative connections 
names = {'left IPG' 'left SMA' 'left dlPFC' 'right IPG', 'right Precuneus', 'right dlPFC'}; % repeating the above, so will have to reorder the nodes again

BMA.Ep(BMA.Pp<.95) = 0;
pos = reshape(BMA.Ep(37:72), 6, 6);
pos = full(pos);
pos = round(pos, 2);
pos(pos>0) = 0;

A = ones(6);
A(pos == 0) = 0;
G = digraph(A, names);

% reordering everything, for intutive graph
order = {'right Precuneus' 'right dlPFC' 'left dlPFC'  'left SMA' 'left IPG' 'right IPG' }
G = reordernodes(G,order);


pos = digraph(pos, names);
pos = reordernodes(pos,order);

LWidths = 5*pos.Edges.Weight/max(pos.Edges.Weight);

EL = pos.Edges.Weight ; % this will be replaced with model number

figure(2)
p = plot(G,'Layout','circle', ...
    'NodeFontSize', 15, nodeColor, 'k', ...
    'EdgeLabel',EL, 'EdgeLabelColor', 'g', 'EdgeFontSize', 12, ...
    'ArrowSize', 10, 'ArrowPosition', 0.99)
