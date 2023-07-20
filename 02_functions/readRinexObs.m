function [errNum,obsOut] = iggReadRinexObs(fileRinex,header,nSat)
% ======================================================================= %
% Function to observations from a RINEX-file
%
%   based on: 
%       'RINEX: The Receiver Independent Exchange Format Version 2.11'
% ----------------------------------------------------------------------- %
% Input:
%   fileRinex ... file identifier from 'iggReadRinex.m'         int[1x1]
%   header ... RINEX header from 'iggReadRinexHeader.m'         struct[1x1]
%   nSat ... number of visible GNSS satellites                  int[1x1]
%
% Output: 
%   errNum ... error code                                       int[1x1]
%   obsOut ... GNSS observations                                double[nxm]
%       n -> number of satellites
%       m -> number of observations taken from RINEX header (sorted in
%            relation to header information)
% ----------------------------------------------------------------------- %
%   @author(s): Florian Zimmermann
%   @last modified: 10.11.2016
%   @contact: zimmermann@igg.uni-bonn.de
% ======================================================================= %

% ----- Initialization of output structure and error code --------------- %
% init error code
errNum = 0;
% output observations structure initialization
obsOut = zeros(nSat,header.nObsTypes);
flagsOut = zeros(nSat,header.nObsTypes);

% ----- Check function input -------------------------------------------- %
if nargin~=3
    errNum = 1005;
    return
end

% ----- Preparation of reading ------------------------------------------ %
% get number maximum number of observation types
nLinesToRead = ceil(header.nObsTypes/5);  % maximum of 5 obs per line
nObsToRead = nLinesToRead * 5;     

% preallocate the line to read
lin = char(32*uint8(ones(16*nObsToRead,1))');

% mask to filter all the possible observations (max 15)
maskObs = false(16,nObsToRead);
maskObs(2:14,:) = true;
maskFlags = true(16,nObsToRead);
maskFlags(2:15,:) = false;
% preallocate a matrix of 15 strings (of length 14 characters)
% notice that each observation element has a max length of 13 char,
% the first character is added as a padding to separate the strings
strObs = char(ones(14,nObsToRead)*32);

% ----- Read observations ----------------------------------------------- %
for i = 1:nSat
    
    % clear line -> fill with spaces
    lin = char(lin*0+32); 
    
    % read all the lines containing the observations (max. 80 characters)
    for l = 1 : (nLinesToRead)
        linTmp = fgetl(fileRinex);
        linLengthTmp = length(linTmp);
        lin((80*(l-1))+(1:linLengthTmp)) = linTmp;  
    end
    
    % convert the lines to a single matrix
    strObs(1:13,:) = (reshape(lin(maskObs(:)),13,nObsToRead));
    strFlags = (reshape(lin(maskFlags(:)),2,nObsToRead));
    
    % read all observations in the string
    fltObs = sscanf(strObs, '%f'); 
    fltFlags = sscanf(strFlags, '%f');
    
    % parsing the observation string
    obsIdx = zeros(header.nObsTypes,1);
    flagIdx = zeros(header.nObsTypes,1);
    for k = 1 : header.nObsTypes
        
        % check if the element is empty
        if (~strcmp(strObs(:,k)','              ')) 
            obsIdx(k) = k;
        end
        if (~strcmp(strFlags(:,k)','  '))
            flagIdx(k) = k;
        end
    end
    
    % assign observations to output structure
    obsOut(i,obsIdx(obsIdx~=0)') = fltObs';
    flagsOut(i,flagIdx(flagIdx~=0)') = fltFlags';
    
end

% check for mapped SNR values and assign them to output structure
if header.snrMapped==1
    [row,~] = find(flagsOut(:));
    tmp = flagsOut(flagsOut~=0);
    dBHz = [25;26.5;29.5;33.5;37;40;43;46.5;49];
    snr = dBHz(tmp);
    obsOut(row+nSat) = snr;
end

return