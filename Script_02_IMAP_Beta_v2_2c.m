% Example 2
% Script #2
%
% Script_02_IMAP_Beta_v2_2c.m
%
% This script reads measured values froma a spreadsheet previously created
% by script #1 (Script_01_IMAP_Beta_v2_2c.m), and then it
% generates violin plots.
%
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
clc

%#########################################################################
%#########################    Edit Starts    #############################
%#########################################################################

% The spreadsheet where you saved the measured values using script #1
bxlsfile  = '/Users/Javier/Desktop/IMAP_output_DATA/IMAP_MEASURED_VALUES_ALLSUBJ_738967673568.xlsx';

% Y-scale limits for violin plots
ylim2plot = [-20 20];

%#########################################################################
%#########################    Edit Ends      #############################
%#########################################################################

% add Violinplot to Matlab's path
IMAP_dir   = fileparts(which(mfilename('fullpath')));
addpath([IMAP_dir filesep 'Violinplot']) ;

% #########################################################################
%          Enter information for plotting using the command window
% #########################################################################

%
% Latencies
%
fprintf('\nMenu: Latencies for ERP mean values :\n\n')
disp('[1] From 300 to 500 ms ')
disp('[2] From 400 to 600 ms ')
disp('[3] From 400 to 800 ms ')
disp('[4] From 400 to 1000 ms ')
disp('[5] From 600 to 800 ms ')
disp('[6] From 600 to 1000 ms ')
% enter a number from the menu
fprintf('\n');
prompt = 'Choose the time window from which you want to plot the mean values (e.g. 3) : ';
xwin = input(prompt);
if isempty(xwin)
      disp('Canceled')
      return
end
if isnumeric(xwin)
      if xwin<1 ||xwin>6 || mod(abs(xwin),1)~=0 || ~isfinite(xwin) || ~isreal(xwin)
            disp('Invalid input. Value must be a positive integer between 1 and 6')
            return
      end
else
      disp('Invalid input. Value must be a positive integer between 1 and 6')
      return
end
switch xwin
      case 1
            latstring = "300500";
      case 2
            latstring = "400600";
      case 3
            latstring = "400800";
      case 4
            latstring = "4001000";
      case 5
            latstring = "600800";
      case 6
            latstring = "6001000";
      otherwise
            error('Invalid latency window (see menu above)')
end

%
% Couples
%
fprintf('\n');
prompt = 'Enter the amount of channel-bin couples (e.g. 3): ';
ncouple = input(prompt);
if isempty(ncouple)
      disp('Canceled')
      return
end
if isnumeric(ncouple)
      if ncouple<1 || mod(abs(ncouple),1)~=0 || ~isfinite(ncouple) || ~isreal(ncouple)
            disp('Invalid input. Value must be a positive integer')
            return
      end
else
      disp('Invalid input. Value must be a positive integer')
      return
end
ways =  [119 104  97 116  32  97 114 101 32 121 111 117  32 115 109 111 107 ...
      105 110 103  63];
if ncouple>=20
      fprintf('\n%s\n',native2unicode(ways));
else
      fprintf('\n');
end
k = 1; stopscript = false;
ch_bin_couple2plot = NaN(ncouple,2);
while k<=ncouple
      ichk = true;
      prompt = sprintf('Enter channel-bin couple #%g  (e.g.  13 3 ) : ', k);
      x = input(prompt,'s');
      x = str2num(x);
      if isnumeric(x)
            if size(x,1)~=1 || size(x,2)~=2
                  ichk = false;
            end
            if ichk
                  if x(1)<=0 || x(2)<=0 || mod(abs(x(1)),1)~=0 || mod(abs(x(2)),1)~=0 ||...
                              ~isfinite(x(1)) || ~isfinite(x(2)) || ~isreal(x(1)) || ~isreal(x(2))
                        ichk = false;
                  end
            end
      else
            ichk = false;
      end
      if ichk
            ch_bin_couple2plot(k,:) = x;
            k = k + 1;
      else
            fprintf('\nInvalid input. Input must be a couple of positive integers\n')
            prompt = 'Do you want to continue? Y/N [Y]: ';
            str = input(prompt,'s');
            if isempty(str)
                  str = 'Y';
            end
            if strcmpi(str,'n')
                  stopscript = true;
                  break
            end
      end
end
if stopscript
      disp('Canceled')
      return
end

% #########################################################################
% #########################################################################
%               Create violin plots (from measured values only)
%
%               ROI: region of interest ;  BOI: bin of interest
% #########################################################################
% #########################################################################

sheets    = sheetnames(bxlsfile);

if ismember(latstring, sheets)
      fprintf('\n\nAll good. I found the sheet "%s" having the mean values for your latency window\n\n', latstring);
else
      fprintf('\n\nToo bad... I could not find the sheet "%s" for your latency window\n\n', latstring);
      disp('Canceled')
      return
end
fprintf('\nI am loading the spreadsheet %s into memory...\n', bxlsfile);
fprintf('Be patient...\n');

T0 = readtable(bxlsfile, 'Sheet', latstring);
T1 = standardizeMissing(T0, '');
T1 = rmmissing(T1);
T2 = convertvars(T1,{...
      'dataset',...
      'channel',...
      'chanLabel',...
      },'categorical');

Tblabel = readtable(bxlsfile, 'Sheet', "BinLabels");
Tblabel = convertvars(Tblabel,{...
      'bini',...
      'binlabel',...
      },'categorical');

fprintf('\n\nI''m now working on the violin plots...one second...\n');

% -------------------------------------------------------------------
% reorganize values into couples ROI&BOI, 400-600 msec
% -------------------------------------------------------------------
labelcouple = cell(1,ncouple);
roiIndx = zeros(1,ncouple);
for cp=1:ncouple
      meanvalues_couples(:,cp)  = T2.(sprintf('Bin%g', ch_bin_couple2plot(cp,2)))(T2.channel==string(ch_bin_couple2plot(cp,1)));
      labelcouple{cp} = sprintf('Ch_%s_%s', string(ch_bin_couple2plot(cp,1)), char(Tblabel.binlabel(ch_bin_couple2plot(cp,2))));
      chIndx(cp)     = ch_bin_couple2plot(cp,1);
end
chIndx = unique(chIndx);

%
% Plot ROI&BOI couples
%
figure
lat1 = latstring{1}(1:3);
lat2 = latstring{1}(4:end);
labelcouple = strrep(labelcouple,'_',' ');
vp1 = violinplot(meanvalues_couples, labelcouple, 'QuartileStyle','shadow'); % call violinplot function
title(sprintf('Violin plot channels [%s]\n', num2str(chIndx)),'interpreter','latex');% ***//EDIT THIS\\***
ylabel(sprintf('Mean values from %s to %s ms ',lat1, lat2));% ***//EDIT THIS\\***
xlim([0, ncouple+1]);
ylim(ylim2plot)

% ---------------------------------------------
disp('Salud!')

