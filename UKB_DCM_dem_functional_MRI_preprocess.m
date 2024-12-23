% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Almgren et al. https://github.com/halmgren/Pipeline_preprint_variability_reliability_DCM_rsfMRI
% This code pre-processes the UKB rs-fmri data
% No slice timing correction needed due to short TR (multiband sequence)

fMRI_datadir = 'C:\Users\Dhiyanesh Srinivasan\OneDrive\Desktop\DATA 2\fMRI';
sMRI_datadir = 'C:\Users\Dhiyanesh Srinivasan\OneDrive\Desktop\DATA 2\sMRI';
subjects = 1;
maxFD = 2.4;

preprocess_UKB_DCM_dem_functional_MRI(subjects, fMRI_datadir, sMRI_datadir, maxFD);

function preprocess_UKB_DCM_dem_functional_MRI(subjects, fMRI_datadir, sMRI_datadir, maxFD)

struct_ID = 'sub-01_T1w.nii'; % updated for your data
funct_ID = 'sub-01_task-rest_bold.nii'; % updated for your data

spm('Defaults', 'fMRI');

for isj = 1:numel(subjects)
    tic;
    EID = subjects(isj);
    matlabbatch = [];

    fprintf(['Running functional preprocessing for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n']);

    % Add this subject's data to the search path
    struct_data_path = sMRI_datadir;
    funct_data_path = fMRI_datadir;
    addpath(struct_data_path);
    addpath(funct_data_path);

    % Verify file existence and print debug information
    fprintf('Checking file paths...\n');
    fprintf('Structural file: %s\n', fullfile(struct_data_path, struct_ID));
    fprintf('Functional file: %s\n', fullfile(funct_data_path, funct_ID));
    if ~exist(fullfile(funct_data_path, funct_ID), 'file')
        error('Functional file not found: %s', fullfile(funct_data_path, funct_ID));
    end
    if ~exist(fullfile(struct_data_path, struct_ID), 'file')
        error('Structural file not found: %s', fullfile(struct_data_path, struct_ID));
    end

    % Skip if subject already pre-processed
    if exist(fullfile(funct_data_path, 'swr_sub-01_task-rest_bold.nii'), 'file') > 0
        fprintf('Done already - skipping...\n');
        continue;
    end

    % Skip if subject partially pre-processed and already aborted due to high motion
    regressors_directory = fullfile(funct_data_path, 'regressors');
    if exist(fullfile(regressors_directory, 'Framewise_Displacement.mat'), 'file')
        load(fullfile(regressors_directory, 'Framewise_Displacement.mat'));
        if max(FD) > maxFD
            clear FD;
            continue;
        end
    end

    % Label data files
    rfMRI_file = fullfile(funct_data_path, funct_ID);
    functional_reference_scan = [rfMRI_file ',1'];
    data = {[functional_reference_scan; (cellstr(spm_select('expand', rfMRI_file)))]};

    structural_reference_scan = fullfile(struct_data_path, [struct_ID ',1']);

    %%%%%%%%%%%%%%%%%%%%%%%
    % Spatial realignment
    %%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{1}.spm.spatial.realign.estwrite.data = data;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r_';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % Save realignment parameters
    rp = load(fullfile(funct_data_path, 'rp_sub-01_task-rest_bold.txt'));
    R = rp;

    mkdir(regressors_directory);
    save(fullfile(regressors_directory, 'rp.mat'), 'R');

    % Calculate and save framewise displacement
    radius = 50; 
    rp(:,4:6) = rp(:,4:6) * radius; 
    rp_diff = diff(rp);
    % Add zero as first element (see Power et al., 2014)
    rp_zero = [zeros(1, size(rp_diff, 2)); rp_diff];
    FD = sum(abs(rp_zero), 2);
    save(fullfile(regressors_directory, 'Framewise_Displacement.mat'), 'FD');

    if max(FD) >= maxFD
        fprintf('Framewise displacement is too high, aborting this subject...\n');
        continue; % Skip this subject if FD too high
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % Coregistration (to first anatomical scan)
    %%%%%%%%%%%%%%%%%%%%%%%
    r_rfMRI_file = fullfile(funct_data_path, ['r_' funct_ID]);
    realigned_volumes = cellstr(spm_select('expand', r_rfMRI_file));

    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {structural_reference_scan};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {functional_reference_scan};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = realigned_volumes;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Segmentation to generate deformation fields
    %%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {structural_reference_scan};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\Program Files\MATLAB\R2024a\spm12\tpm\TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Normalize
    %%%%%%%%%%%%%%%%%%%%%%%
    deformation_field = {fullfile(struct_data_path, ['y_' struct_ID])};
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = deformation_field;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = realigned_volumes;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];

    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2.4, 2.4, 2.4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Smooth
    %%%%%%%%%%%%%%%%%%%%%%%
    wr_rfMRI_file = fullfile(funct_data_path, ['wr_' funct_ID]);
    normalised_volumes = cellstr(spm_select('expand', wr_rfMRI_file));

    matlabbatch{1}.spm.spatial.smooth.data = normalised_volumes;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    rmpath(funct_data_path); % To conserve memory
    rmpath(struct_data_path); % To conserve memory

    fprintf(['COMPLETED functional preprocessing for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n']);
    fprintf(['Run time: ' num2str(toc) ' seconds\n']);
end

end

function extract_timeseries(subjects, fMRI_datadir, maxFD)
    funct_ID = 'sub-01_task-rest_bold.nii'; % updated for your data

    % ROI labels and MNI coordinates
    [ROI_list] = specify_ROI();

    spm('Defaults', 'fMRI');

    failed_subjs = [];
    retained_subjs = [];
    for isj = 1:numel(subjects)
        tic
        EID = subjects(isj);

        fprintf(['Extracting ROI timeseries for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n']);

        funct_data_path = fullfile(fMRI_datadir, num2str(EID));
        regressors_directory = fullfile(funct_data_path, 'regressors');
        SPM_file_directory = fullfile(funct_data_path, 'SPM');  % Ensure this is correctly defined

        % Ensure the SPM directory exists
        if ~exist(SPM_file_directory, 'dir')
            mkdir(SPM_file_directory);
        end

        % Skip if timeseries extraction is already complete for this subject
        if exist(fullfile(SPM_file_directory, 'T_ext_flag.mat'), 'file')
            load(fullfile(SPM_file_directory, 'T_ext_flag.mat'));
            if strcmp(T_ext_flag, 'failure')
                failed_subjs = [failed_subjs; EID];
            else
                retained_subjs = [retained_subjs; EID];
            end
            clear T_ext_flag
            fprintf('Skipping subject - already done\n');
            continue
        end

        % Skip subject if motion was too excessive
        fprintf('Loading Framewise_Displacement.mat...\n');
        load(fullfile(regressors_directory, 'Framewise_Displacement.mat'), 'FD');
        if max(FD) >= maxFD
            fprintf('Skipping subject - motion too excessive\n');
            failed_subjs = [failed_subjs; EID];
            continue
        end

        % Try/catch statement to pick up subjects where preprocessing/timeseries extraction failed
        try 
            fprintf('Starting timeseries extraction process...\n');
            matlabbatch = [];
            % Add this subject's pre-processed fMRI data to the MATLAB search path
            addpath(funct_data_path);
            pre_proc_rfMRI = fullfile(funct_data_path, 'swr_sub-01_task-rest_bold.nii');

            % Extract TR from subject header. Should be 0.735 for UKB
            fprintf('Extracting TR from NIfTI header using spm_vol...\n');
            V = spm_vol(pre_proc_rfMRI);
            TR = V(1).private.timing.tspace;

            if round(TR, 3) ~= 0.735
                warning('TR read from header is %s...UKB TR should be 0.735', num2str(round(TR, 3)));
            end

            % Create regressors (DCT, motion, CSF, WM)
            fprintf('Creating regressors...\n');
            [nifti_images, DCT] = UKB_DCM_BuildRegressors(funct_data_path, pre_proc_rfMRI, TR);

            fprintf('Creating SPM directory and files...\n');
            % Create directory
            if ~exist(SPM_file_directory, 'dir')
                mkdir(SPM_file_directory);
            end

            matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_file_directory};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
                fullfile(regressors_directory, 'dct_set.mat'), % Discrete cosine basis set
                fullfile(regressors_directory, 'VOI_CSF_signal_1.mat'), % CSF signal
                fullfile(regressors_directory, 'VOI_WM_signal_1.mat'), % WM signal
                fullfile(regressors_directory, 'rp.mat') % 6 realignment (motion) parameters
            };

            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

            %% fMRI model estimation
            fprintf('Starting fMRI model estimation...\n');
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            %% Define contrasts
            fprintf('Defining contrasts...\n');
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Effects_of_interest';
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(size(DCT, 2)); % F contrast over whole cosine set
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.delete = 0;

            spm_jobman('run', matlabbatch);

            % Extract Timeseries from each ROI
            thresholds = nan(1, size(ROI_list, 1));

            for VOI_number = 1:length(ROI_list) % Loop over ROIs
                clear matlabbatch;
                fprintf('Extracting timeseries for ROI %d...\n', VOI_number);
                cd(SPM_file_directory);
                matlabbatch{1}.spm.util.voi.spmmat(1) = {['SPM.mat']};
                matlabbatch{1}.spm.util.voi.adjust = 1;
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.name = ROI_list{VOI_number, 1};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat(1) = {['SPM.mat']};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = ROI_list{VOI_number, 2};
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10; % Creates 10 mm sphere centered around pre-defined coordinates
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = ROI_list{VOI_number, 2};
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8; % Creates 8 mm sphere centered around peak of F-contrast within the previous 10 mm sphere
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';

                % Voxels only included if they exceed the threshold (i1) AND lie
                % within the fixed 10 mm sphere (i2) AND lie within the mobile 8 mm
                % sphere (i3)
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';

                count = 0;
                err_count = 0;
                while count == err_count % Make threshold more liberal if no surviving voxels (up to max of 0.05)
                    fprintf('Attempt %d for ROI %d with threshold %f...\n', count+1, VOI_number, matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh);
                    if count == 0
                        matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
                    elseif count == 1
                        cd(SPM_file_directory);
                        matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01;
                    elseif count == 2
                        cd(SPM_file_directory);
                        matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
                    elseif count >= 3
                        matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = nan; % These sessions are excluded
                        break;
                    end
                    spm_jobman('run', matlabbatch);

                    if ~exist(fullfile(SPM_file_directory, ['VOI_' ROI_list{VOI_number, 1} '_1.mat']), 'file')
                        fprintf('No surviving voxels for ROI %d at threshold %f\n', VOI_number, matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh);
                        err_count = err_count + 1;
                    end
                    count = count + 1;
                end

                thresholds(VOI_number) = matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
            end

            save('thresholds.mat', 'thresholds');

            clear matlabbatch;

            rmpath(funct_data_path); % To conserve memory

            fprintf(['COMPLETED timeseries extraction for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n']);
            fprintf(['Run time: ' num2str(toc) ' seconds\n']);

            % Save flag to show that timeseries extraction is complete
            T_ext_flag = 'success'; 
            save(fullfile(SPM_file_directory, 'T_ext_flag.mat'), 'T_ext_flag');
        catch ME
            fprintf(['FAILED timeseries extraction for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n']);
            fprintf('Error message: %s\n', ME.message);
            fprintf(['Run time: ' num2str(toc) ' seconds\n']);

            % Save flag to show that timeseries extraction is complete
            T_ext_flag = 'failure';
            save(fullfile(SPM_file_directory, 'T_ext_flag.mat'), 'T_ext_flag');
        end
    end
end

function [ROI_list] = specify_ROI()
    % Define ROI labels and MNI coordinates
    ROI_list = {
        'CSF_signal', [0, 0, 0];
        'WM_signal', [0, 0, 0];
        % Add more ROIs as needed
    };
end
