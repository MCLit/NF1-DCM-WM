% code for creating VOIs - 6mm spheres located at centres of group clusters
clear
clc

root = '/Users/user/data';
spmpath= '/Users/user/Documents/MATLAB/spm12/';

cd(root);
spm('Defaults','fMRI');

subcode = [];% list of subject folder codes
subcode = num2str(subcode, '%05.f');

voi_name = 'gg_left_dlPFC';
voi_coordinates = [-42 20 25];
inner = 6; % size of the VOI radius

for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']); % location of SPM file with first level results
    matlabbatch{1}.spm.util.voi.adjust = 1; % this regresses out the effects of motion etc. read more at:
    % 35.3.8 Extracting data from regions of https://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf
    % and Adjustment in https://jacoblee.net/occamseraser/2018/01/03/extracting-rois-for-ppi-analysis-using-spm-batch/index.html
    matlabbatch{1}.spm.util.voi.session = 1; % run 1 
    matlabbatch{1}.spm.util.voi.name = voi_name; % what to name the volume
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates; % location of the center of the sphere
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner; % radius of the sphere
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1; % whether sphere is fixed
    matlabbatch{1}.spm.util.voi.expression = 'i1'; % logical and/or statements
    spm_jobman('run',matlabbatch);
    clear job
end

%%
clc
voi_name = 'gg_right_dlPFC';
voi_coordinates = [42 32 30];
for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']);

    matlabbatch{1}.spm.util.voi.adjust = 1; 
    matlabbatch{1}.spm.util.voi.session = 1; % run 1
    matlabbatch{1}.spm.util.voi.name = voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run',matlabbatch);
    clear job
end


%%
clc
voi_name = 'gg_left_IPG';
voi_coordinates = [-28 -50 38];

for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']);

    matlabbatch{1}.spm.util.voi.adjust = 1;
    matlabbatch{1}.spm.util.voi.session = 1; % run 1
    matlabbatch{1}.spm.util.voi.name = voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run',matlabbatch);
    clear job
end

%%
clc
voi_name = 'gg_right_IPG';
voi_coordinates = [30 -58 42];

for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']);
    matlabbatch{1}.spm.util.voi.adjust = 1; 
    matlabbatch{1}.spm.util.voi.session = 1; % run 1
    matlabbatch{1}.spm.util.voi.name = voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run',matlabbatch);
    clear job
end

%%
clc
voi_name = 'gg_left_SMA';
voi_coordinates = [-8 2 62];

for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']);

    matlabbatch{1}.spm.util.voi.adjust = 1;
    matlabbatch{1}.spm.util.voi.session = 1; % run 1
    matlabbatch{1}.spm.util.voi.name = voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run',matlabbatch);
    clear job
end

%%
clc
voi_name = 'gg_right_Precuneus';
voi_coordinates = [8 -62 48];

for s = 1:size(subcode)
    cd([root,filesep,subcode(s,:)]);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr([root,filesep,subcode(s,:),filesep,'results/SPM.mat']);

    matlabbatch{1}.spm.util.voi.adjust = 1;
    matlabbatch{1}.spm.util.voi.session = 1; % run 1
    matlabbatch{1}.spm.util.voi.name = voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = voi_coordinates;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = inner;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
    spm_jobman('run',matlabbatch);
    clear job
end


