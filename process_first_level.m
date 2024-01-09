clc
clear

root = '/Users/user/data';
spmpath= '/Users/user/Documents/MATLAB/spm12/';

cd(root);
spm('Defaults','fMRI');

folders = dir('0*');% you can change this if the files are named in another way
for i=1:length(folders)
    tmp(i,:)=folders(i,1).name;
end
subcode = (tmp);% you can manually remove cases that do not have good data


for s = 1:size(subcode,1)

    cd([root,filesep,subcode(s,:)]);
    mkdir 'results'
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr([root,filesep,subcode(s,:), filesep, 'results']);
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.5;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    f  = dir([root subcode(s,:) '/Run1_combined/swu*.img']);
    for i = 1:144
        all{i,:} = [f(i).folder '/' f(i).name];
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = all;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = '0-back';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = [2.038780846
        62.03278182
        122.0411404
        182.0311636
        242.0363961
        302.0411573];
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 24;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = '2-back';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = [32.04023511
        92.03502516
        152.0243816
        212.0414791
        272.0216088
        332.0266615];
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 24;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});



    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr([root,filesep,subcode(s,:), filesep, 'Run1_combined', filesep, 'art_regression_outliers_and_movement_u0001.mat']);
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';


matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;




matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'n-back';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 0
                                                        0 1];
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '0-back';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = '2-back';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = '0>2-back';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = '2>0-back';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
    spm_jobman('run',matlabbatch);
    disp('first level complete')
    disp(s)
    cd(root)

end