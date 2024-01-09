%%
clc
clear
root = '/Users/user/data';

cd(root);
spm('Defaults','fMRI');

folders = dir('0*');% you can change this if the files are named in another way
for i=1:length(folders)
    tmp(i,:)=folders(i,1).name;
end
subcode = (tmp);% you can manually remove cases that do not have good data

n_scan=144;%input('enter number of dymanics per run'); % number of dynamic scans
n_run=1;%input('enter number of runs');% number of runs
slicetime='Y';%input('Do you want to do slice time correction (e.g. suitable for event-related designs)?: Y/N?','s');


parrec_n1='Run1';%input('Enter 1st par rec filename ','s');
allruns=parrec_n1;%enter par rec file name
do_sub=1:length(subcode);

%% cleaning - BE CAREFUL WHAT YOU CLEAN
% for s = 1:size(subcode)
%     cd(root);
% 
%     cd(subcode(s,:));
% 
% 
%     delete('T1/*')
%     delete('Run1_combined/*')
%     delete('Run1_short/*')
%     delete('Run1_long/*')
% 
%     delete('results/*')
%         cd ..
% end

%%
for s = 1:size(subcode)
    cd(subcode(s,:));
    % extract T1 into nii
    % mkdir('T1');
    % 
    % dicm2nii('T1.par',pwd,0);
     spm_parrec2nifti([pwd, '/T1/T1.PAR']) % best method
    % 
    % nii_hdr2nii([pwd, '/T1.hdr']) % dicm2nii not working, just need to make nii out of hdr

    movefile('T1*','T1'); %spm wants nii file in subject folder for
cd(root)
end
% NOW: Check Reg - is your orgin at anterior commissure?

%%
for s = 1:size(subcode)
    cd(root)
    cd(subcode(s,:));

    % extract PARREC to short and long echo, takes into account if you have
    % multiple runs - you need to name each run at beginning.
    pname_a = [root,filesep,subcode(s,:),filesep];
    for i=1:n_run
        fname_a = [allruns(i,:),'.REC'];
        dual_echo_analyse(pname_a,fname_a);
    end

    % realignment
    clear job

    tmpcwdfolder = [root,filesep,subcode(s,:),filesep,parrec_n1,'_short',filesep];
    tmpcwdfolder2 = [root,filesep,subcode(s,:),filesep,parrec_n1,'_long',filesep];

    % Do slice timing correction if requested
    if slicetime=='Y'
        AHslicetimecorrection_nofieldmap(n_run,root,subcode(s,:),parrec_n1);
    end

    % select img that are slice corrected
    if n_run==1
        run1s = cellstr(spm_select('FPList', tmpcwdfolder, '^a0.*\.img$'));
        run1l = cellstr(spm_select('FPList', tmpcwdfolder2, '^a0.*\.img$'));
    end

    % realign fMRI dynamics
    if n_run==1
        AHrealign_orig_multruns(n_run,root,subcode(s,:),parrec_n1);
        % AHrealign_orig_multruns(n_run,root,subcode,parrec_n1);
        % AHrealign_orig_multruns(n_run,run1s,run1l);
    end

    % combine short and long echo images, for all runs
    cd([root,filesep,subcode(s,:),filesep]);
    for i=1:n_run
        mkdir([allruns(i,:),'_combined']);
        tmpcwdfolder = [root,filesep,subcode(s,:),filesep,allruns(i,:),'_short',filesep];
        tmpcwdfolder2 = [root,filesep,subcode(s,:),filesep,allruns(i,:),'_long',filesep];
        if slicetime=='Y'
            urun1 = cellstr(spm_select('FPList', tmpcwdfolder, '^ra.*\.img$'));
            urun2 = cellstr(spm_select('FPList', tmpcwdfolder2, '^ra.*\.img$'));
        end
        for r=1:length(urun1) % for each scan
            %load long and short volumes and turn into vector
            vol1=spm_vol(urun1{r});
            vol1=reshape(spm_read_vols(vol1),1,[]);
            vol2=spm_vol(urun2{r});
            vol2=reshape(spm_read_vols(vol2),1,[]);
            %average vector
            comb1=(vol1+vol2)./2;
            %create combined volume
            tmp_img=spm_vol(urun1{r});
            tmp_img.fname=[root,filesep,subcode(s,:),filesep,allruns(i,:),'_combined',filesep,sprintf('u%s.img',num2str(r,'%04i'))];
            tmp_img.pinfo=[0;0;0];
            spm_write_vol(tmp_img,reshape(comb1,tmp_img.dim(1),tmp_img.dim(2),tmp_img.dim(3)));
        end
        %need a mean image too for DPARFA
        mean1 = cellstr(spm_select('FPList', tmpcwdfolder, '^meana.*\.img$'));
        mean2 = cellstr(spm_select('FPList', tmpcwdfolder2, '^meana.*\.img$')); % changed to be coreg version for long
        vol1=spm_vol(mean1{1});
        vol1=reshape(spm_read_vols(vol1),1,[]);
        vol2=spm_vol(mean2{1});
        vol2=reshape(spm_read_vols(vol2),1,[]);
        comb1=(vol1+vol2)./2;
        tmp_img=spm_vol(mean1{1});
        tmp_img.fname=[root,filesep,subcode(s,:),filesep,allruns(i,:),'_combined',filesep,sprintf('meana%s.img',num2str(1,'%04i'))];
        tmp_img.pinfo=[0;0;0];
        spm_write_vol(tmp_img,reshape(comb1,tmp_img.dim(1),tmp_img.dim(2),tmp_img.dim(3)));
    end

    clear job

end
%% coregistration, check outputs prior and post to ensure that coreg is working fine
for s = 1:size(subcode)
    % coregister
    tmpT1= [root,filesep,subcode(s,:),filesep,'T1/T1.nii,1'];
    cc=1;
    % take all other functional images with Coreg fucntion
    if n_run==1
        tmpcwdfolder = [root,filesep,subcode(s,:),filesep,parrec_n1,'_combined',filesep];
        run1 = cellstr(spm_select('FPList', tmpcwdfolder, '^u.*\.img$'));
        job{cc}.spm.spatial.coreg.estimate.other = {run1{:}, [root,filesep,subcode(s,:),filesep,parrec_n1,'_long',filesep,'meana0001.img,1'],[root,filesep,subcode(s,:),filesep,parrec_n1,'_combined',filesep,'meana0001.img,1']}'; %changed to coreg long mean as well as need to create combined mean for DPARSFA
    end

    job{cc}.spm.spatial.coreg.estimate.ref = {tmpT1}; % Specfiy T1, because of tilted brain, we try to adjust source to reference space
    job{cc}.spm.spatial.coreg.estimate.source = {[root,filesep,subcode(s,:),filesep,parrec_n1,'_short',filesep,'\meana0001.img,1']};
    %     job{cc}.spm.spatial.coreg.estimate.other = {}; %specified above when accounting for run
    job{cc}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    job{cc}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    job{cc}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    job{cc}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    %run job list
    spm_jobman('run',job);
    disp('co-registration done')
    disp(s)
    cd ..

end

%% segmentation
for s = 1:size(subcode)


    cd ([root subcode(s,:) '/T1/'])

    clear job

    T1_img = ({[root subcode(s,:) '/T1/T1.nii,1']});
    job{1}.spm.spatial.preproc.channel.vols = T1_img;
    job{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    job{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    job{1}.spm.spatial.preproc.channel.write = [0 0];
    job{1}.spm.spatial.preproc.tissue(1).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,1'};
    job{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    job{1}.spm.spatial.preproc.tissue(1).native = [1 1];
    job{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    job{1}.spm.spatial.preproc.tissue(2).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,2'};
    job{1}.spm.spatial.preproc.tissue(2).ngaus =1;
    job{1}.spm.spatial.preproc.tissue(2).native = [1 1];
    job{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    job{1}.spm.spatial.preproc.tissue(3).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,3'};
    job{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    job{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    job{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    job{1}.spm.spatial.preproc.tissue(4).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,4'};
    job{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    job{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    job{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    job{1}.spm.spatial.preproc.tissue(5).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,5'};
    job{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    job{1}.spm.spatial.preproc.tissue(5).native = [0 0];
    job{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    job{1}.spm.spatial.preproc.tissue(6).tpm = {'/Users/user/Documents/MATLAB/spm12/tpm/TPM.nii,6'};
    job{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    job{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    job{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    job{1}.spm.spatial.preproc.warp.mrf = 1;
    job{1}.spm.spatial.preproc.warp.cleanup = 1;
    job{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    job{1}.spm.spatial.preproc.warp.affreg = 'mni';
    job{1}.spm.spatial.preproc.warp.fwhm = 0;
    job{1}.spm.spatial.preproc.warp.samp = 3;
    job{1}.spm.spatial.preproc.warp.write = [0 0];
    job{1}.spm.spatial.preproc.warp.vox = NaN;
    job{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];

    spm_jobman('run',job);
    disp('segmentation done')
    disp(s)
    cd(root)
    
    
end
%% DARTEL - create templates that optimise WM and GM alignment

for s = 1:size(subcode)
    cd ([root subcode(s,:) '/T1'])

    clear job

    rc1_img = ({[root subcode(s,:) '/T1/rc1T1.nii,1']});
    rc2_img = ({[root subcode(s,:) '/T1/rc2T1.nii,1']});

    job{1, 1}.spm.tools.dartel.warp.images{1,1} = [rc1_img; rc2_img];
    job{1, 1}.spm.tools.dartel.warp.settings.template = 'template';
    job{1, 1}.spm.tools.dartel.warp.settings.rform = 0;
    job{1, 1}.spm.tools.dartel.warp.settings.param(1).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(1).rparam = [4,2,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(1).K = 0;
    job{1, 1}.spm.tools.dartel.warp.settings.param(1).slam = 16;


    job{1, 1}.spm.tools.dartel.warp.settings.param(2).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(2).rparam = [2,1,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(2).K = 0;
    job{1, 1}.spm.tools.dartel.warp.settings.param(2).slam = 8;


    job{1, 1}.spm.tools.dartel.warp.settings.param(3).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(3).rparam = [1,0.500000000000000,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(3).K = 1;
    job{1, 1}.spm.tools.dartel.warp.settings.param(3).slam = 4;


    job{1, 1}.spm.tools.dartel.warp.settings.param(4).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.500000000000000,0.250000000000000,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(4).K = 2;
    job{1, 1}.spm.tools.dartel.warp.settings.param(4).slam = 2;


    job{1, 1}.spm.tools.dartel.warp.settings.param(5).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.250000000000000,0.125000000000000,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(5).K = 4;
    job{1, 1}.spm.tools.dartel.warp.settings.param(5).slam = 1;


    job{1, 1}.spm.tools.dartel.warp.settings.param(6).its = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.250000000000000,0.125000000000000,1.000000000000000e-06];
    job{1, 1}.spm.tools.dartel.warp.settings.param(6).K = 6;
    job{1, 1}.spm.tools.dartel.warp.settings.param(6).slam = 0.500000000000000;


    job{1, 1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.0100;
    job{1, 1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
    job{1, 1}.spm.tools.dartel.warp.settings.optim.its = 3;


    spm_jobman('run',job);
    disp('DARTEL template ready')
    disp(s)
    cd(root)
end

%% DARTEL - normalise to MNI (note: specific for fMRI!)

for s = 1:size(subcode)
    cd ([root subcode(s,:) '/T1'])

    clear job

    tmplt = ({[root subcode(s,:) '/T1/template_6.nii']});
    job{1, 1}.spm.tools.dartel.mni_norm.template = tmplt; 


    flws1 = ({[root subcode(s,:) '/T1/u_rc1T1_template.nii']});
    job{1, 1}.spm.tools.dartel.mni_norm.data.subj.flowfield = flws1;

    f  = dir([root subcode(s,:) '/Run1_combined/u*.img']);
    for i = 1:144
        all{i,:} = [f(i).folder '/' f(i).name];
    end
    job{1, 1}.spm.tools.dartel.mni_norm.data.subj.images = all;
    job{1, 1}.spm.tools.dartel.mni_norm.vox = [2.5 2.5 2.5]; % we can change this to be up or downsampled
    job{1, 1}.spm.tools.dartel.mni_norm.bb = [NaN, NaN, NaN; NaN, NaN, NaN];
    job{1, 1}.spm.tools.dartel.mni_norm.preserve = 0;
    job{1, 1}.spm.tools.dartel.mni_norm.fwhm = [6,6,6];

    spm_jobman('run',job);
    disp('DARTEL normalisation & smoothing done')
    disp(s)
    cd(root)
end