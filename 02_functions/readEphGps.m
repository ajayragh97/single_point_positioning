function [errNum, ephGps, ephTime, ephIdx] = iggReadEphGps(filename,projectname,maxSat)
% ======================================================================= %
% Function to read the Broadcast ephemeris for GPS satellites
%
%   based on:
%       'RINEX: The Receiver Independent Exchange Format Version 2.11'
% ----------------------------------------------------------------------- %
% Input:
%   filename ... Ephemeris filename                            str[-]
%   projectname ...  project foldername                        str[-]
%   maxSat ... maximum number of satellites per GNSS           int[1x1]
%
% Output:
%   errNum ... error code                                      int[1x1]
%   ephGps ... GPS ephemeris                               double(35xn)
%       n -> number of epochs with ephemeris updates
%   ephTime ... ephemeris update times [SOW]                   double[kxmaxSat]
%   ephIdx ... columns in ephGps related to each satellite     int[kxmaxSat]
%       k -> max number of updates +1
% ----------------------------------------------------------------------- %
%   @author(s): F. Zimmermann, C. Eling
%   @last modified: 23.03.2017
%   @contact: zimmermann@igg.uni-bonn.de
% ======================================================================= %

% ----- Initialization of output structure and error code --------------- %
% init error return
errNum = 0;

% initialize output structure
ephGps = [];
ephTime = [];
ephIdx = [];

% ----- Check function input -------------------------------------------- %
if nargin~=3
    errNum = 1007;
    return
end

if ~ischar(filename) || ~ischar(projectname)
    errNum = 1008;
    return
end

% ----- Move through file to look for end of header --------------------- %
fid = fopen(filename, 'r');
iter = 1;
while 1
    % get next line of file
    tline = fgetl(fid);
    % find end of header
    if (strfind(tline,'END OF HEADER'))
        line_start = iter;
        break
    end
    iter = iter+1;
end
fclose(fid);

% ----- Move through file to count number of epochs to read ------------- %
% initialize epoch vector with enough space for possible entries
line_sat = zeros(500,1);

fid = fopen(filename, 'r');
j = 1;
iteration = 1;
dummy = textscan(fid,'%s %s %s',1,'headerlines',line_start-1);
tline = fgetl(fid);
while 1
    % get next line of file
    tline = fgetl(fid);
    % find end of file
    if ~ischar(tline)
        break
    end
    % get line number of new epoch
    if ~strcmp(tline(1:3),'   ')
        line_sat(j) = iteration+line_start;
        j = j+1;
    end
    % increment counter
    iteration = iteration+1;
end
fclose(fid);

% delete unused entries in epoch vector
line_sat(line_sat==0) = [];

% ----- Read ephemeries data -------------------------------------------- %
% preallocate output structure
ephGps = zeros(35,length(line_sat));

for k = 1:length(line_sat)
    % open file
    fid = fopen(filename,'r');
    % read data for satellite k
    prn = textscan(fid,'%3c',1,'headerlines',line_sat(k)-1);
    prnStr = regexp(prn{1},'\d+','match');
    data = textscan(fid,'%n',34);
    % close file
    fclose(fid);
    % store data in ephemeris matrix
    ephGps(:,k)=[str2double(prnStr{1});data{1}];
end

% ----- Build ephemeris time and time-index matrices -------------------- %
% get maximal number of ephemeris updates for one satellite
mUp = max(histc(ephGps(1,:),1:maxSat));
% preallocate time and index matrix
ephTime = zeros(mUp+1,maxSat);
ephIdx = zeros(size(ephTime));

% fill time and index matrix
for i = 1:maxSat
    % get all update times for satellite i
    satIdx = ephGps(1,:) == i;
    satTime = ephGps(2:7,satIdx);
    satTime(1,:) = satTime(1,:)+2000;
    
    % convert times to SOW
    sow = zeros(size(satTime,2),1);
    for j = 1:size(satTime,2)
        sow(j) = getSOW(satTime(:,j));
    end
    
    % fill time matrix
    ephTime(2:length(sow)+1,i) = sow;
    % fill index matrix
    ephIdx(2:length(sow)+1,i) = find(satIdx');
end

return