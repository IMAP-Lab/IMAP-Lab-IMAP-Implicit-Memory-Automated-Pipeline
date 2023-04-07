% Author: Javier Lopez-Calderon
% Newencode Analytics
% www.newencode.com
% Talca, Chile
% 2022-2023
%
% This function gets all the filenames within a folder (pathname) with an
% specific extension (dataExt). Additionally, it searches for a common string
% among the found filenames located at the end of their names. 
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

function [Subjects, commontail] = getFileNamesEEG(pathname, dataExt)

dsetList = dir(fullfile(pathname, ['*' dataExt]));
Subjects = {dsetList.name}';
Subjects = strrep(Subjects, dataExt,''); % remove extension (temporary)
nsubj    = length(Subjects);
Lmax     = max(cellfun(@length, Subjects));
grandMat = zeros(nsubj-1, Lmax);

if nsubj>1
        SubjStr = upper(string(Subjects));
        for j=2:nsubj
                StrTemp0 = char(SubjStr(1));
                N0 = length(StrTemp0);
                StrTempK = char(SubjStr(j));
                Nk = length(StrTempK);
                if Nk<Lmax
                        StrTempK = [repmat(char(randi([97 122])),1,Lmax-Nk) StrTempK];
                        Nk = length(StrTempK);
                end

                if N0<Lmax
                        StrTemp0 = [repmat(char(randi([97 122])),1,Lmax-N0) StrTemp0];
                end
                Q = StrTempK( mod(bsxfun(@plus,(0:Nk-1)',0:Nk-1),Nk)+1 );
                logMaskStr = bsxfun(@eq, StrTemp0, Q);
                [~, indxrow] = max( sum(logMaskStr, 2) );
                grandMat(j-1,:) = logMaskStr(indxrow, :);
        end

        fm    = logical(prod(grandMat,1));
        nfm   = length(fm);
        tmask = false(1,nfm);
        k=0;
        while all(fm(end-k:end))
                tmask(1,nfm-k) = true;
                k=k+1;
        end
        commontail = StrTemp0(tmask);
else
        commontail = [];
end

