% Example 1
% Script #1
%
% Script_01_IMAP_Beta_v2_2c.m
%
% This script does the following (*):
%
% - Applies the function RA_equaltrials_Beta_GH001.m to your
%   previously epoched datasets.
% - Creates both individual (equated)averaged ERPs and the grand averaged ERP.
% - Exports counting stats per both confidence and status levels.
% - Also, all measured values are exported to Excel
%
% * This script is based on data/paradigm/codes of Addante,(2015) Neuroimage
% % and Addante, Lopez-Calderon et al (2023)
% ----------------------------------------------------------------------
%                            Authors:
% Javier Lopez-Calderon, PhD              Richard J. Addante, PhD
% Newencode Analytics                     Florida Institute of Technology
% Talca, Chile                            Melbourne, FL
% ----------------------------------------------------------------------
% Created: April-September, 2022
% Updated (help) for GitHub: December 2022-January 2023
% ----------------------------------------------------------------------
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% ----------------------------------------------------------------------
% www.newencode.com
% 2023
% ----------------------------------------------------------------------

clear
close all
clc

%#########################################################################
%#########################    Edit Starts    #############################
%#########################################################################

% Note for User: You will need to change each of the items 1:5 listed below
% in lines 48:58 (e.g. your own file folders for data input & output), and
% possibly the event codes listed in point #6.

% 1. Datasets folder path. Put in this folder all of the *.set EEG files you want to be processed by IMAP
pathname_read_EEG        = '/Users/Javier/Desktop/TEST_pre_Github';

% 2. Output EEG folder path. IMAP will save here the processed EEG data
pathname_output_proc_EEG = fullfile('/Users/Javier/Desktop/TEST_pre_Github/EQUALTRIALS_EPOCH/IMAP_output_EEG');

% 3. Output ERP folder path. IMAP will save here the newly created EEG data
pathname_output_proc_ERP = fullfile('/Users/Javier/Desktop/TEST_pre_Github/EQUALTRIALS_EPOCH/IMAP_output_ERP');

% 4. Output report folder (texts, grand averaged ERP, spreadsheets, etc)
% This is the folder that will be used by script #2 in order to read the
% generated spreadsheet
path_report = '/Users/Javier/Desktop/TEST_pre_Github/IMAP_output_DATA';

% 5. Add EEGLAB and requiered plugins (optional)
addEeglab = false; % true or false
% Eeglab's full path (optional)
eeglab_dir = '/Users/Javier/Dropbox (Personal)/eeglab2023.0';

% 6. Next Section That User May Need to Update: Your event codes for bins
% Numeric event codes, found in epoched datasets, corresponding to each
% confidence level (e.g. 1 to 5) and their statuses (e.g. old vs new)

% if your task has a larger or smaller number of confidence levels, simply
% copy (increasing its index) or delete the correponding (two) lines per
% level specified above (use numeric event codes triggered by your own task)

% event codes for confidence level 1, statuses old and new
ConfidenceLevelCode(1).old = [6111 6112 6113 6114 6115 6511 6512 6513 6514 6515]; %bin1 = old1
ConfidenceLevelCode(1).new = [6311 6312 6313 6314 6315]; %bin2 = new1

% event codes for confidence level 2, statuses old and new
ConfidenceLevelCode(2).old = [6121 6122 6123 6124 6125 6521 6522 6523 6524 6525];%bin3 = old2
ConfidenceLevelCode(2).new = [6321 6322 6323 6324 6325];%bin4 = new2

% event codes for confidence level 3, statuses old and new
ConfidenceLevelCode(3).old = [6131 6132 6133 6134 6135 6531 6532 6533 6534 6535];%bin5 = old3
ConfidenceLevelCode(3).new = [6331 6332 6333 6334 6335];%bin6 = new3

% event codes for confidence level 4, statuses old and new
ConfidenceLevelCode(4).old = [6141 6142 6143 6144 6145 6541 6542 6543 6544 6545];%bin7 = old4
ConfidenceLevelCode(4).new = [6341 6342 6343 6344 6345];%bin8 = new4

% event codes for confidence level 5, statuses old and new
ConfidenceLevelCode(5).old = [6151 6152 6153 6154 6155 6551 6552 6553 6554 6555];%bin9 = old5
ConfidenceLevelCode(5).new = [6351 6352 6353 6354 6355];%bin10 = new5

% minimum number of trials accepted (for equated number of trials)
MinNumTrial = 12;

% log file name (to be saved)
fileOutStat = 'IMAP_Trial_Counting_ALLSUBJ.xlsx';

% measured values file name (to be saved)
fileOutMea = 'IMAP_MEASURED_VALUES_ALLSUBJ.xlsx';

% Optional: low-pass filter for individual ERP data (Saved separately)
erp_lpf_cutoff = 30;  % cutoff in Hz (e.g. 30). No filter -> 0

% Optional: low-pass filter for grand averaged ERP
% (Saved separately. For instance, for paper's figures)
gavg_lpf_cutoff = 30;  % cutoff in Hz (e.g. 30). No filter -> 0

% IMPORTANT!(if you don't include a channel or bin here then it won't be
% plotted or exported for stats). Check your channel and bin definitions
% at "chanFormulas" and "binFormulas" below.
binArray   = 1:21; % bin(s) index to be included for exporting values (and plotting later)
chanArray  = 1:43; % chan(s)index to be included for exporting values (and plotting later)

%#########################################################################
%#########################    Edit Ends      #############################
%#########################################################################

% eventCodeArray     - It's a 1xN cell array of event codes, where N is the number of
%                    "Confidence" responses (e.g. 5 responses).
%                    * Each element is a new 1xM cell array where M is
%                      the different "Status" we are comparing accross within
%                      a given "Confidence" response (in this case, they represent
%                      "old" versus "new").
%                    * Each element of this cell array is an array that contains the event
%                      codes that belong to each "Status".
nconflevel = length(ConfidenceLevelCode);
eventCodeArray = cell(1, nconflevel);

for jc=1:nconflevel
        eventCodeArray{1,jc} = {ConfidenceLevelCode(jc).old  ConfidenceLevelCode(jc).new};
end
if exist(pathname_output_proc_EEG, 'dir')==0
        % if folder doesn't exist then make it
        mkdir(pathname_output_proc_EEG);
end
if exist(pathname_output_proc_ERP, 'dir')==0
        % if folder doesn't exist then make it
        mkdir(pathname_output_proc_ERP);
end
if exist(path_report, 'dir')==0
        % if folder doesn't exist then make it
        mkdir(path_report);
end

% -------------------------------------------------------------------------
% add EEGLAB to Matlab's path
if addEeglab

        addpath(eeglab_dir, genpath([eeglab_dir filesep 'functions'])) ;
        rmpath(genpath([eeglab_dir filesep 'functions' filesep 'octavefunc'])) ;

        % add plugins folder
        pdir = dir([eeglab_dir filesep 'plugins']) ;
        pdirfindex = [pdir.isdir];
        pdir = {pdir.name};
        plugin_dir = pdir(pdirfindex & startsWith(string(pdir), lettersPattern));
        plugin_dir = strcat(eeglab_dir, filesep, 'plugins', filesep, plugin_dir, ';') ;
        addpath([plugin_dir{:}]) ;

        % add ERPLAB to Matlab's path
        if exist('erplab_default_values', 'file')
                erplab_dir = which('eegplugin_erplab.m') ;
                erplab_dir = erplab_dir(1:strfind(erplab_dir, 'eegplugin_erplab.m')-1) ;
                addpath(genpath(erplab_dir)) ;
        else; error('Please make sure ERPLAB is on your path') ;
        end
end

% -------------------------------------------------------------------------
% PREPARE DATASET LIST TO BE PROCESSED-------------------------------------
% IMAP only needs the data folder path to automatically load the datasets to be included
% in the current analysis
dataExt  = '.set';
[Subjects, commontail] = getFileNamesEEG(pathname_read_EEG, dataExt);
nsubj    = length(Subjects);

% -------------------------------------------------------------------------
ERPSubjectList = cell(10,1); % whole path and filename of good subjects
pointGoodErp = 1;
pointBadErp  = 1;

%
% Main loop
%
for k=1:nsubj

        %
        % filename to be loaded (without extension .set)
        %
        sname = Subjects{k};

        fprintf('\n---------------------------------------------\n');
        fprintf('Loading %s.set\n', sname);

        %         try
        %
        % Load dataset (epoch dataset. Bad epochs must be already marked at this point)
        %
        EEG = pop_loadset('filename',[sname '.set'],'filepath', pathname_read_EEG);

        %
        % apply RA_equaltrials_Beta001GH.m
        %
        EEG = RA_equaltrials_Beta001GH(EEG, eventCodeArray);

        %
        % Save new dataset (replace common "tail" in filenames if any)
        %
        if ~isempty(commontail)
                sname = strrep(sname, commontail, '_IMAP_OUT_EEG');
        else
                sname = [sname '_IMAP_OUT_EEG'];
        end
        fprintf('\n###Saving dataset...\n');
        EEG.setname = sname;
        EEG = pop_saveset( EEG,  'filename', [EEG.setname  '.set'], 'filepath', pathname_output_proc_EEG);

        %
        % Get "good" epoch indices per bin to attach this info to the
        % ERP structure
        goodepochsAll = bitand(not(sum(EEG.reject.rejmanualE,1)>0), not(EEG.reject.rejmanual));
        nbin = EEG.EVENTLIST.nbin;
        FinalSelectedEpochs = zeros(nbin, EEG.trials);
        for ibin = 1:nbin
                FinalSelectedEpochs(ibin, epoch4bin(EEG, ibin)) = 1;
                FinalSelectedEpochs(ibin, :) = bitand(FinalSelectedEpochs(ibin, :), goodepochsAll);
        end

        %
        % Get averaged ERPs
        %
        ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
        % attach good epoch indices per bin
        ERP.ntrials.FinalSelectedEpochs = FinalSelectedEpochs;

        %
        % "chanFormulas" : Add new channels (ROIs)
        %
        chanFormulas = {...
                'ch31 = (ch13+ch12+ch18)/3 label = Posterior3_P3_Pz_P4',...
                'ch32 = (ch12+ch13+ch18+ch11+ch21+ch15+ch17)/7 label = Posterior7_Pz_P3_P4_Cp1_Cp2_O1_O2',...
                'ch33 = (ch12+ch13+ch18+ch11+ch21)/5 label = Posterior5_Pz_P3_P4_Cp1_Cp2',...
                'ch34 = (ch12+ch13+ch18+ch15+ch16+ch17)/6 label = Posterior6_Pz_P3_P4_O1_O2_Oz',...
                'ch35 = (ch14+ch13+ch12+ch18+ch19+ch15+ch16+ch17+ch10+ch11+ch20+ch21)/12 label = Posterior11_P7_P3_Pz_P4_P8_O1_Oz_O2_Cp5_Cp1_Cp2_Cp6',...
                'ch36 = (ch12+ch13+ch18+ch15+ch16+ch17+ch14+ch19)/8 label = Posterior8_P7_P3_Pz_P4_P8_O1_Oz_O2',...
                'ch37 = (ch1+ch2+ch3+ch4+ch6+ch7+ch26+ch27+ch28+ch29+ch30)/11 label = Frontal11_Fp1_Fp2_F7_F8_F3_F4_Fz_Fc5_Fc1_Fc2_Fc6',...
                'ch38 = (ch3+ch4+ch28+ch29+ch2+ch6+ch7+ch26+ch27)/9 label = Frontal8_F7_F8_F3_F4_Fz_Fc5_Fc1_Fc2_Fc6',...
                'ch39 = (ch2+ch7+ch27+ch3+ch28)/5 label = Frontal5Mid_Fz_Fc1_Fc2_F3_F4',...
                'ch40 = (ch2+ch7+ch27)/3 label = Frontal3_Fz_Fc1_Fc2',...
                'ch41 = (ch2+ch3+ch28)/3 label = Frontal3_Fz_F3_F4',...
                'ch42 = (ch18+ch21+ch12)/3 label = PosteriorLeft_P4_CP2_Pz',...
                'ch43 = (ch13+ch11+ch12)/3 label = PosteriorRight_P3_CP1_Pz',...
                'ch44 = (ch12+ch21+ch11)/3 label = PosteriorCentral_Pz_Cp1_Cp2'...
                };

        ERP = pop_erpchanoperator( ERP, chanFormulas, 'KeepLocations',  1, 'Warning', 'off' );

        %
        %  "binFormulas" : Add new bins
        %
        binFormulas = {...
                'bin11 = wavgbin(1,3,5,7,9) label IMAP Old12345',...
                'bin12 = wavgbin(2,4,6,8,10) label IMAP New12345',...
                'bin13 = bin11-bin12 label IMAP Old-New',...
                'bin14 = wavgbin(1,3,5,7) label IMAP Old1234',...
                'bin15 = wavgbin(2,4,6,8) label IMAP New1234',...
                'bin16 = bin14-bin15 label IMAP Old1234-New1234',...
                'bin17 = wavgbin(1,3,5) label IMAP Old123',...
                'bin18 = wavgbin(2,4,6) label IMAP New123',...
                'bin19 = bin17-bin18 label IMAP Old123-New123',...
                'bin20 = bin14-bin12 label IMAP Old1234-New12345',...
                'bin21 = bin17-bin12 label IMAP Old123-New12345',...
                'bin22 = bin11-bin18 label IMAP 12345-equated123',...
                'bin23 = wavgbin(1,3,5,9) label Old1235'...
                'bin24 = wavgbin(2,4,6,10) label New1235'...
                'bin25 = bin23-bin24 label Old1235vsNew1235'...
                'bin26 = bin23-bin12 label Old1235vsNew12345'...
                };

        ERP = pop_binoperator( ERP, binFormulas);

        %
        % Save erpset
        %
        fprintf('\n###Saving erpset...\n');
        fprintf('\nSaving subject %s \n', sname);
        erpname = sprintf('%s_ERP', sname); % subject's filename
        ERP     = pop_savemyerp(ERP, 'erpname', erpname, 'filename', [erpname '.erp'], 'filepath', pathname_output_proc_ERP);
        fprintf('Processed erpset was saved!\n\n');

        if all(ERP.ntrials.accepted(11:12) >= MinNumTrial)
                fprintf('\n\n%s\n\n', repmat('$',1,100));
                fprintf('*** YAY! Subject %s was included in the current study \n\n', sname);
                fprintf('%s\n', repmat('$',1,100));
                % whole path and filename of good subjects
                ERPSubjectList{pointGoodErp} = fullfile(pathname_output_proc_ERP, [erpname '.erp']);

                %
                % Store ERP in ALLERP
                %
                ALLERP(pointGoodErp) = ERP;

                %
                % Store amount of trials per confidence and status
                %
                ARTDETECT_ACCE(pointGoodErp).erpname = ERP.erpname;
                for hh=1:length([ERP.ntrials.accepted])
                        condname = replace(ERP.bindescr{hh}, {' ','-','+','*','/'},'_');
                        ARTDETECT_ACCE(pointGoodErp).(condname)  = ERP.ntrials.accepted(hh);
                end
                pointGoodErp = pointGoodErp + 1;
        else
                fprintf('\n\n%s\n\n', repmat('@',1,100));
                fprintf('*** Boo... Subject %s was excluded from the current study... Sorry. \n\n', sname);
                fprintf('%s\n', repmat('@',1,100));
                %
                % Store amount of trials per confidence and status
                %
                ARTDETECT_EXCLU(pointBadErp).erpname = ERP.erpname;
                for hh=1:length([ERP.ntrials.accepted])
                        condname = replace(ERP.bindescr{hh}, {' ','-','+','*','/'},'_');
                        ARTDETECT_EXCLU(pointBadErp).(condname)  = ERP.ntrials.accepted(hh);
                end
                pointBadErp = pointBadErp + 1;
        end

        % *************************************************************
        % EEG data filtering (optional. See header for freq cutoff
        % values)
        % *************************************************************
        if erp_lpf_cutoff>0
                % LOW-PASS FILTER erp_lpf_cutoff Hz (ERPLAB)-----------------------------------
                fprintf('Lowpass filtering the data...\n');
                ERP = pop_filterp( ERP,  1:ERP.nchan , 'Cutoff', erp_lpf_cutoff, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  2 );

                %
                % Save filtered erpset
                %
                fprintf('\n###Saving filtered erpset version...\n');
                erpname = sprintf('%s_LOWPASS%g', erpname, erp_lpf_cutoff);
                ERP     = pop_savemyerp(ERP, 'erpname', erpname, 'filename', [erpname '.erp'], 'filepath', pathname_output_proc_ERP);
                fprintf('Filtered erpset was saved!\n\n');
        end
        clear EEG ERP
        %         catch
        %                 fprintf('\n\n\n***Oops! %s crashes...\n\n', sname);
        %         end
end

fprintf('\n\n*** Still working...Please wait...\n');

%
% Export subjects list to text
%
Terp = cell2table(ERPSubjectList);
slistfile = fullfile(path_report, 'ERP_IMAP_SubjectList.txt');
writetable(Terp, slistfile, 'WriteVariableNames', false);

%
% Get Grand Average
%
if pointGoodErp>1
        ERP = pop_gaverager( ALLERP , 'Criterion',50, 'Erpsets', 1:length(ALLERP), 'ExcludeNullBin', 'on', 'SEM', 'on' );

        chanLabels = {ERP.chanlocs.labels}';
        binLabels  = ERP.bindescr';

        %
        % Save Grand Average
        %
        fprintf('*** Saving grand averaged ERP...\n');
        erpname = 'GAVG_ERP_IMAP';
        ERP  = pop_savemyerp(ERP, 'erpname', erpname, 'filename', [erpname '.erp'], 'filepath', path_report);

        if gavg_lpf_cutoff>0
                LPFALLERP = ALLERP;
                for jj=1:length(ALLERP)
                        LPFALLERP(jj) = pop_filterp( ALLERP(jj),  1:ERP.nchan , 'Cutoff', gavg_lpf_cutoff, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  2 );
                end
                LPFGAVG = pop_gaverager( LPFALLERP , 'Criterion',50, 'Erpsets', 1:length(LPFALLERP), 'ExcludeNullBin', 'on', 'SEM', 'on' );
                %
                % Save low-pass filtered Grand Average
                %
                fprintf('*** Saving lowpass filtered grand averaged ERP...\n');
                lpgavgerpname = sprintf('GAVG_ERP_IMAP_%g_Hz_LowPass', gavg_lpf_cutoff);
                pop_savemyerp(LPFGAVG, 'erpname', lpgavgerpname, 'filename', [lpgavgerpname '.erp'], 'filepath', path_report);
        end

else
        % script 1 won't continue if no good ERP remains
        fprintf('\n\n\n*** Grand Average failed ***...\n\n');
        return
end

%
% Export trial stats from confidence and status conditions
%
if pointGoodErp>1
        TA = struct2table(ARTDETECT_ACCE)
        writetable(TA,fullfile(path_report, ['Good_' fileOutStat]));
end

if pointBadErp>1
        TE = struct2table(ARTDETECT_EXCLU)
        writetable(TE,fullfile(path_report, ['Exclu_' fileOutStat]));
end

% #########################################################################
%                              Extract measured values
% #########################################################################

%
% Extract and export to text mean values from 400-600 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_400_600.txt');
[~, meanvalues400600] = pop_geterpvalues( slistfile, [ 400 600],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

%
% Extract and export to text mean values from 300-500 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_300_500.txt');
[~, meanvalues300500] = pop_geterpvalues( slistfile, [ 300 500],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

%
% Extract and export to text mean values from 400-800 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_400_800.txt');
[~, meanvalues400800] = pop_geterpvalues( slistfile, [ 400 800],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

%
% Extract and export to text mean values from 600-1000 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_600_1000.txt');
[~, meanvalues6001000] = pop_geterpvalues( slistfile, [ 600 1000],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

%
% Extract and export to text mean values from 400-1000 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_400_1000.txt');
[~, meanvalues4001000] = pop_geterpvalues( slistfile, [ 400 1000],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

%
% Extract and export to text mean values from 600-800 (bin x chan x subj)
%
meafile = fullfile(path_report, 'MeanValues_600_800.txt');
[~, meanvalues600800] = pop_geterpvalues( slistfile, [ 600 800],  binArray,  chanArray , 'Baseline', 'pre',...
        'Binlabel', 'on', 'FileFormat', 'wide', 'Filename', meafile, 'InterpFactor',  1,...
        'Measure', 'meanbl', 'Mlabel', 'mean',  'Resolution',  3 );

% #########################################################################
%                           Export measured values
%
%           There will be one spreadsheet (*.xlsx) with multiple sheets
%
% #########################################################################

%
% reorganize values into a 3D array (sub x bin x chan)
%
meanvalues400600_reorg  = permute(meanvalues400600, [3 1 2]);
meanvalues300500_reorg  = permute(meanvalues300500, [3 1 2]);
meanvalues400800_reorg  = permute(meanvalues400800, [3 1 2]);
meanvalues6001000_reorg = permute(meanvalues6001000,[3 1 2]);
meanvalues4001000_reorg = permute(meanvalues4001000,[3 1 2]);
meanvalues600800_reorg  = permute(meanvalues600800, [3 1 2]);

%
% Output spreadsheet final name
%
taildate  = sprintf('%12.0f',now*10^6); % add date and time suffix to create unique filename
[~,fileOutMea] = fileparts(fileOutMea);
tablename = fullfile(path_report, sprintf('%s_%s.xlsx',fileOutMea, taildate));

% ----------------------------------------------------------------------------
% reorganize 400600 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m400600 = permute(meanvalues400600_reorg,[1 3 2]);
m400600 = reshape(m400600,[],size(meanvalues400600_reorg,2),1);
chindx  = repmat(chanArray',size(meanvalues400600_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues400600_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues400600_reorg,3),1),length(chindx),1);

%
% create 400600 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm400600 = array2table(m400600,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm400600.('dataset') = datasetname;
Tm400600.('channel') = chindx;
Tm400600.('chanLabel') = chanLab;
Tm400600 = movevars(Tm400600,{'dataset','channel', 'chanLabel'},'Before',1);

% create 400600 spreadsheet
writetable(Tm400600, tablename ,'Sheet',sprintf('400600'),'Range','A1');

% ----------------------------------------------------------------------------
% reorganize 300500 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m300500 = permute(meanvalues300500_reorg,[1 3 2]);
m300500 = reshape(m300500,[],size(meanvalues300500_reorg,2),1);
chindx  = repmat(chanArray',size(meanvalues300500_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues300500_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues300500_reorg,3),1),length(chindx),1);

%
% create 300500 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm300500 = array2table(m300500,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm300500.('dataset') = datasetname;
Tm300500.('channel') = chindx;
Tm300500.('chanLabel') = chanLab;
Tm300500 = movevars(Tm300500,{'dataset','channel', 'chanLabel'},'Before',1);

% create 300500 spreadsheet
writetable(Tm300500, tablename ,'Sheet',sprintf('300500'),'Range','A1');

% ----------------------------------------------------------------------------
% reorganize 400800 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m400800 = permute(meanvalues400800_reorg,[1 3 2]);
m400800 = reshape(m400800,[],size(meanvalues400800_reorg,2),1);
chindx  = repmat(chanArray',size(meanvalues400800_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues400800_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues400800_reorg,3),1),length(chindx),1);

%
% create 400800 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm400800 = array2table(m400800,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm400800.('dataset') = datasetname;
Tm400800.('channel') = chindx;
Tm400800.('chanLabel') = chanLab;
Tm400800 = movevars(Tm400800,{'dataset','channel', 'chanLabel'},'Before',1);

% create 400800 spreadsheet
writetable(Tm400800, tablename ,'Sheet',sprintf('400800'),'Range','A1');

% ----------------------------------------------------------------------------
% reorganize 6001000 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m6001000 = permute(meanvalues6001000_reorg,[1 3 2]);
m6001000 = reshape(m6001000,[],size(meanvalues6001000_reorg,2),1);
chindx   = repmat(chanArray',size(meanvalues6001000_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues6001000_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues6001000_reorg,3),1),length(chindx),1);

%
% create 6001000 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm6001000 = array2table(m6001000,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm6001000.('dataset') = datasetname;
Tm6001000.('channel') = chindx;
Tm6001000.('chanLabel') = chanLab;
Tm6001000 = movevars(Tm6001000,{'dataset','channel', 'chanLabel'},'Before',1);

% create 6001000 spreadsheet
writetable(Tm6001000, tablename ,'Sheet',sprintf('6001000'),'Range','A1');

% ----------------------------------------------------------------------------
% reorganize 4001000 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m4001000 = permute(meanvalues4001000_reorg,[1 3 2]);
m4001000 = reshape(m4001000,[],size(meanvalues4001000_reorg,2),1);
chindx  = repmat(chanArray',size(meanvalues4001000_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues4001000_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues4001000_reorg,3),1),length(chindx),1);

%
% create 4001000 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm4001000 = array2table(m4001000,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm4001000.('dataset') = datasetname;
Tm4001000.('channel') = chindx;
Tm4001000.('chanLabel') = chanLab;
Tm4001000 = movevars(Tm4001000,{'dataset','channel', 'chanLabel'},'Before',1);

% create 4001000 spreadsheet
writetable(Tm4001000, tablename ,'Sheet',sprintf('4001000'),'Range','A1');

% ----------------------------------------------------------------------------
% reorganize 600800 values into a 3-column, 2D array (colums: sub, bin, chan)
% ----------------------------------------------------------------------------
m600800 = permute(meanvalues600800_reorg,[1 3 2]);
m600800 = reshape(m600800,[],size(meanvalues600800_reorg,2),1);
chindx  = repmat(chanArray',size(meanvalues600800_reorg,1),1);
chanLab  = repmat(chanLabels(chanArray),size(meanvalues600800_reorg,1),1);
datasetname = reshape(repmat(string(ERP.workfiles),size(meanvalues600800_reorg,3),1),length(chindx),1);

%
% create 600800 Table and then save it into a specific sheet within a
% single spreadsheet ("fileOutMea" variable at header)
%
Tm600800 = array2table(m600800,'VariableNames', append("Bin", strtrim(string(num2str([binArray]'))')));
Tm600800.('dataset') = datasetname;
Tm600800.('channel') = chindx;
Tm600800.('chanLabel') = chanLab;
Tm600800 = movevars(Tm600800,{'dataset','channel', 'chanLabel'},'Before',1);

% create 600800 spreadsheet
writetable(Tm600800, tablename ,'Sheet',sprintf('600800'),'Range','A1');

%##########################################################################
% write bin labels on the same spreadsheet
Binx.bini = binArray';
Binx.binlabel = {binLabels{binArray}}';
writetable(struct2table(Binx), tablename ,'Sheet',"BinLabels",'Range','A1');

% ---------------------------------------------
disp('Salud!')
