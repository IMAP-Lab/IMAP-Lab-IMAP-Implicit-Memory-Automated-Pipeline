% SYNTAX:
% EEG = RA_equaltrials_Beta001GH(EEG, eventCodeArray, minNumEp);
%
% INPUT:
%
% EEG             - bin-based epoched dataset (artifact detection/correction already done)
% eventCodeArray  - It's a 1xN cell array of event codes, where N is the number of
%                   "Confidence" responses (e.g. 5 responses).
%                    * Each element is a new 1xM cell array where M is
%                      the different "Status" we are comparing accross within
%                      a given "Confidence" response (in this case, they represent 
%                      "old" versus "new").
%                    * Each element of this cell array is an array that contains the event 
%                      codes that belong to each "Status".
%                    * See the following example:
%
% eventCodeArray = {...
%     {[6111 6112 6113 6114 6115 6511 6512 6513 6514 6515] [6311 6312 6313 6314 6315]}...
%     {[6121 6122 6123 6124 6125 6521 6522 6523 6524 6525] [6321 6322 6323 6324 6325]}...
%     {[6131 6132 6133 6134 6135 6531 6532 6533 6534 6535] [6331 6332 6333 6334 6335]}...
%     {[6141 6142 6143 6144 6145 6541 6542 6543 6544 6545] [6341 6342 6343 6344 6345]}...
%     {[6151 6152 6153 6154 6155 6551 6552 6553 6554 6555] [6351 6352 6353 6354 6355]}...
%     };
%
% >> eventCodeArray =
% 
%   1×5 cell array
% 
%     {1×2 cell}    {1×2 cell}    {1×2 cell}    {1×2 cell}    {1×2 cell}
%
%
% minNumEp        - minimum admissible amount of trials per "Status"
%                   Optional. Default is 1
%
%
% OUTPUT:
%
% EEG             - bin-based epoched dataset (ready for getting averaged
%                   ERPs!)
%
%                   IMPORTANT: 
%                   1) This bin-based epoched dataset has fewer
%                   epochs than the original one (duh! of course!). 
%                   2) It also has new bin labels/descriptions according to
%                   the Confidence and Status condition/response.
%
%                   Let's see an example. We have 5 confidence responses
%                   for 2 statuses ('old' and 'new')
%
%                   Confidence   |  Status  |  bin label |  (new)code label
%                   ------------------------------------------------------
%                        1       |     1    |    Bin 1   |     101
%                        1       |     2    |    Bin 2   |     102
%                        2       |     1    |    Bin 3   |     201
%                        2       |     2    |    Bin 4   |     202
%                        3       |     1    |    Bin 5   |     301
%                        3       |     2    |    Bin 6   |     302
%                        4       |     1    |    Bin 7   |     401
%                        4       |     2    |    Bin 8   |     402
%                        5       |     1    |    Bin 9   |     501
%                        5       |     2    |    Bin 10  |     502
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% Current function works with epoched-data (EEGLAB's struct).
% See explanation here below:
%
%#########################################################################
%         ###################  How it works ######################
%#########################################################################
%
% Relevant EEG's fields for this function:
% 
% EEG.epoch--> EEG.EVENTLIST.eventinfo--> EEG.EVENTLIST.bdf
%
% -length of EEG.epoch = as many as number of epochs after epoching data.
%
% -length of EEG.EVENTLIST.eventinfo = as many as event code in the
%                                      continuous file.
%
% -length of EEG.EVENTLIST.bdf = as many as bins described in the
%                                bin descriptor file (bdf).
%
% *** Example. At Subject 61_Exp019_ICAclean_updated.set, the epoch # 10 has
%         eventitem = 28 (from all events in the continuous file). So,
%
% >> EEG.epoch(10)
%
% ans =
%
%   struct with fields:
%
%              event: 12
%        eventbepoch: {[10]}
%          eventbini: {[457]}
%      eventbinlabel: {'B457(S63)'}
%       eventbvmknum: {[28]}
%        eventbvtime: {[]}
%       eventchannel: {[0]}
%     eventcodelabel: {'S63'}
%      eventduration: {[15.6250]}
%        eventenable: {[1]}
%          eventflag: {[0]}
%          eventitem: {[28]}  <---------- item 28
%       eventlatency: {[0]}
%          eventtype: {'B457(S63)'}
%
%
% Thus, the eventitem # 28 has the following info:
%
% >> EEG.EVENTLIST.eventinfo(28)
%
% ans =
%
%   struct with fields:
%
%          item: 28
%          code: 63
%      binlabel: 'B457(S63)'  <----------- bin 457
%     codelabel: 'S63'
%          time: 84.7720
%        spoint: 2.1703e+04
%          dura: 4
%          flag: 0
%        enable: 1
%          bini: 457          <---------- bin 457
%        bepoch: 10
%       channel: 0
%        bvtime: [0×1 double]
%       bvmknum: 28
%
%  So we can see here that the bin number (index) -that captured this
%  epoch- was the bin #457 ("bini" and "binlabel"),and the description about
%  this bin is at EEG.EVENTLIST.bdf(457)
%
%  >> EEG.EVENTLIST.bdf(457)
%
% ans =
%
%   struct with fields:
%
%      expression: '.{63}{2}{2}'
%     description: '6322'      <---------- description from bdf (4 digits!)
%         prehome: []
%          athome: [1×1 struct]
%        posthome: [1×2 struct]
%         namebin: 'BIN 457'
%          rtname: []
%         rtindex: []
%              rt: []
%
% Therefore, the original 4-digit event code, that identifies the type of stimulus,
% can be found at:
%
% EEG.EVENTLIST.bdf(457).description
% ans =
%      '6322'
%
% ########################################################################
%
% Example of use:
%
% eventCodeArray = {...
%     {[6111 6112 6113 6114 6115 6511 6512 6513 6514 6515] [6311 6312 6313 6314 6315]}...
%     {[6121 6122 6123 6124 6125 6521 6522 6523 6524 6525] [6321 6322 6323 6324 6325]}...
%     {[6131 6132 6133 6134 6135 6531 6532 6533 6534 6535] [6331 6332 6333 6334 6335]}...
%     {[6141 6142 6143 6144 6145 6541 6542 6543 6544 6545] [6341 6342 6343 6344 6345]}...
%     {[6151 6152 6153 6154 6155 6551 6552 6553 6554 6555] [6351 6352 6353 6354 6355]}...
%     };
% 
% EEG = pop_loadset('filename','63_Exp019_ICAclean_updated.set','filepath','C:\Users\raddante\Desktop\ICA Cleaned_Updated');
% EEG = RA_equaltrials_Beta001GH(EEG, eventCodeArray);
%
% ----------------------------------------------------------------------
% Author:
% Javier Lopez-Calderon
% Newencode Analytics
% www.newencode.com
% Talca, Chile
% April-September, 2022
%
% Inspired by Carter Luck's equaltrials.m function
%
% For: Richard Addante's Lab, Florida Tech
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
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function EEGout = RA_equaltrials_Beta001GH(EEGin, eventCodeArray, minNumEp)
if nargin<1
      doc RA_equaltrials_Beta001GH
      help RA_equaltrials_Beta001GH
end
if nargin<2
      error('RA_equaltrials_V2:EVENTLIST','You must specify a cell array of event codes. See help RA_equaltrials_Beta001GH\n')
end
if nargin<3
      minNumEp = 1;
end
if isempty(minNumEp)
      minNumEp = 1;
end
if minNumEp<0
      error('RA_equaltrials_V2:EVENTLIST','"minNumEp" must be a positive integer! (or zero)\n')
end
if ~iscell(eventCodeArray)
      error('RA_equaltrials_V2:EVENTLIST','"eventCodeArray" must be a cell array of integer numbers\n')
end
if ~isfield(EEGin,'EVENTLIST')
      error('RA_equaltrials_V2:EVENTLIST','\nOops! "EVENTLIST" structure was not found. Run ERPLAB''s binlister first\n')
end
%
% initialize  EEGout
%
EEGout = EEGin;
EEGout = rmfield(EEGout,'data');
EEGout = rmfield(EEGout,'epoch');
EEGout = rmfield(EEGout,'event');
if isfield(EEGout,'urevent')
      EEGout = rmfield(EEGout,'urevent');
end
if isfield(EEGout,'EVENTLIST')
      EEGout = rmfield(EEGout,'EVENTLIST');
end
ntrials = EEGin.trials; % number of epochs at EEGin
nConfidence   = length(eventCodeArray);

% initialize  selected epochs
selectedEpochs = NaN(ntrials,3);

% Iterate through all epochs
for EpochIndex = 1:ntrials
      
      nStatus = zeros(1,nConfidence); % initialize amount of Status within current Confidence
      
      % iterate through all Confidences
      for ConfidenceIndex = 1:nConfidence % e.g. 1:5
            
            % get Status values within current Confidence
            StatusEventCode = eventCodeArray{1,ConfidenceIndex};
            
            % update amount of Status within current Confidence
            nStatus(ConfidenceIndex) = length(StatusEventCode);
            
            % find time-locked event ("bin" index) within current epoch (multiple events within an epoch)
            if iscell(EEGin.epoch(EpochIndex).eventlatency)
                    eventbini = cell2mat(EEGin.epoch(EpochIndex).eventbini(ismember([EEGin.epoch(EpochIndex).eventlatency{:}],0)));
            else % single event within an epoch
                    eventbini = EEGin.epoch(EpochIndex).eventbini;
            end

            % get string having 4-digit format. E.g. '6135'
            evencode4digit = str2num([EEGin.EVENTLIST.bdf(eventbini).description]);
            
            % select epoch: add 'old'/'new' and 'Confidence' identity to the current epoch
            for StatusIndex = 1:nStatus(ConfidenceIndex)
                  if any(ismember(evencode4digit ,StatusEventCode{1,StatusIndex}))% belong to old/new?
                        selectedEpochs(EpochIndex,1) = ConfidenceIndex ;  % Confidence
                        selectedEpochs(EpochIndex,2) = StatusIndex ;% e.g. 1=old; 2=new
                        selectedEpochs(EpochIndex,3) = EpochIndex;  % epoch index
                  end
            end
      end
end
% remove NaNs
selectedEpochs(any(isnan(selectedEpochs), 2), :) = [];

% initialize pointers
indxnewepo = 1; % pointer for including selected epoch

for ConfidenceIndex = 1:nConfidence % Confidence loop
      % initialize masks
      countermask       = false(nStatus(ConfidenceIndex), size(selectedEpochs,1));
      counterConfidence = zeros(1,nStatus(ConfidenceIndex));
      for StatusIndex = 1:nStatus(ConfidenceIndex) %Status loop
            % find epoch indices belonging to Confidence "ConfidenceIndex" AND Status "StatusIndex"
            countermask(StatusIndex,:) = selectedEpochs(:,1)==ConfidenceIndex & selectedEpochs(:,2)==StatusIndex;
            % count successfully found epoch indices
            counterConfidence(StatusIndex)   = nnz(countermask(StatusIndex,:));
      end
      %
      % find min number of epochs between/among Status
      %
      [minNumberOfEpochs, minStatus] = min(counterConfidence);
      if minNumberOfEpochs >= minNumEp
            % copy selected epochs into EEGout
            for k=1:nStatus(ConfidenceIndex)
                  % epoch indices belonging to Confidence "ConfidenceIndex" AND Status k
                  epindx = selectedEpochs(countermask(k,:), 3);
                  if ~isempty(epindx)
                        if k~=minStatus
                              randidx = randperm(length(epindx));
                              epindx  = epindx(randidx(1:minNumberOfEpochs));
                        end
                        % selected epochs' loop
                        for h=1:minNumberOfEpochs
                              %
                              % EEGout's event's fields
                              %
                              EEGout.event(indxnewepo).epoch  = indxnewepo;
                              EEGout.event(indxnewepo).item   = indxnewepo;
                              EEGout.event(indxnewepo).bepoch = indxnewepo;
                              if rem(k, 2) == 0 % even
                                    binix  = ConfidenceIndex*k ;% new bini ("bin index")
                              else % odd
                                    binix  = ConfidenceIndex*k + abs(ConfidenceIndex-k);% new bini ("bin index")
                              end
                              EEGout.event(indxnewepo).bini      = binix;
                              EEGout.event(indxnewepo).binlabel  = sprintf('B%g(%g)',binix, ConfidenceIndex*100 + k);
                              EEGout.event(indxnewepo).codelabel = sprintf('%g', ConfidenceIndex*100 + k);
                              EEGout.event(indxnewepo).type      = EEGout.event(indxnewepo).binlabel;
                              EEGout.event(indxnewepo).enable    = 1;
                              EEGout.event(indxnewepo).flag      = 0;
                              [~, sampzero]                      = closest(EEGin.times, 0);
                              EEGout.event(indxnewepo).latency   = sampzero+ EEGin.pnts*(indxnewepo-1); % crucial!                              
                              %
                              % EEGout's epochs' fields
                              %
                              EEGout.epoch(indxnewepo)               = EEGin.epoch(epindx(h));
                              EEGout.epoch(indxnewepo).eventtype     = ConfidenceIndex*100 + k; % 101, 102, 201, etc
                              EEGout.data(:,:,indxnewepo)            =  EEGin.data(:,:,epindx(h));
                              EEGout.epoch(indxnewepo).event         = indxnewepo;
                              EEGout.epoch(indxnewepo).eventbepoch   = {[indxnewepo]};
                              EEGout.epoch(indxnewepo).eventbini     = {[binix]};
                              EEGout.epoch(indxnewepo).eventbinlabel = {EEGout.event(indxnewepo).binlabel};
                              EEGout.epoch(indxnewepo).eventbvmknum  = {[]};
                              EEGout.epoch(indxnewepo).eventbvtime   = {[]};
                              EEGout.epoch(indxnewepo).eventchannel  = {[]};
                              EEGout.epoch(indxnewepo).eventcodelabel= {EEGout.event(indxnewepo).codelabel};
                              EEGout.epoch(indxnewepo).eventduration = {[5]};
                              EEGout.epoch(indxnewepo).eventenable   = {[1]};
                              EEGout.epoch(indxnewepo).eventflag     = {[0]};
                              EEGout.epoch(indxnewepo).eventitem     = {[indxnewepo]};
                              EEGout.epoch(indxnewepo).eventlatency  = {[0]};
                              EEGout.epoch(indxnewepo).eventtype     = ConfidenceIndex*100 + k;                                                                 
                              % EEGout's epochs' pointer
                              indxnewepo = indxnewepo + 1;
                        end
                  else
                        fprintf('WARNING: Status #%d  from  Confidence  #%d  does not have matching event codes...\n',k , ConfidenceIndex);
                  end
            end
      else          
            fprintf('WARNING:  Confidence  #%d  has a minimum number of trials of %d...\n', ConfidenceIndex, minNumberOfEpochs);
      end
end
if indxnewepo==1
      error('RA_equaltrials_V2:EVENTLIST','\nOops! I couldn''t find matching event codes...Please check your dataset\n')
end
% rebuild a "dummy" EVENTLIST
EEGout.EVENTLIST      = EEGin.EVENTLIST;
EEGout.EVENTLIST.trialsperbin = [];
EEGout.EVENTLIST.elname = [];
p=1;
binlabel = cell(1,nConfidence*mode(nStatus));
for a=1:nConfidence
      for b=1:mode(nStatus)
            binlabel{p} = sprintf('B%g(%g)',p, a*100 + b);
            p = p + 1;
      end
end
EEGout.EVENTLIST.nbin = length(binlabel);
for k=1:EEGout.EVENTLIST.nbin
      EEGout.EVENTLIST.bdf(k).expression  =  'NA';
      bindesc = char(regexp(binlabel{k}, '\((\d*)\)', 'tokens', 'once'));
      bindesc = sprintf('Confidence %s Status %s',bindesc(1), bindesc(2:3));
      % modification (workaround)---
      bindesc = strrep(bindesc,'Status 01','Status Old');
      bindesc = strrep(bindesc,'Status 02','Status New');
      % -----------------------------     
      EEGout.EVENTLIST.bdf(k).description = bindesc;
      EEGout.EVENTLIST.bdf(k).prehome     = 'NA';
      EEGout.EVENTLIST.bdf(k).athome      = 'NA';
      EEGout.EVENTLIST.bdf(k).posthome    = 'NA';
      namebinx     =  char(regexp(binlabel{1}, '(\w*)(', 'tokens', 'once'));
      EEGout.EVENTLIST.bdf(k).namebin     = strrep(namebinx,'B','Bin');
end
EEGout.EVENTLIST.eldate = char(datetime('now','TimeZone','local','Format', 'MMM dd, yyyy HH:mm:ss.SSS'));
EEGout.EVENTLIST.eventinfo = EEGout.event;
EEGout.trials = size(EEGout.data, 3);
% EEGout.urevent = EEGin.urevent;
% update EEG.reject.rejmanual
EEGout.reject.rejmanualE = single(zeros(EEGin.nbchan, EEGout.trials));
EEGout.reject.rejmanual  = single(zeros(1 , EEGout.trials));
% check everything
EEGout = eeg_checkset(EEGout);
disp('Done')
