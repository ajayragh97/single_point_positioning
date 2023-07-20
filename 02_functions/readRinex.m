function [errNum, header, obs, date, sow, obsOut] = iggReadRinex(filename,projectname,maxSat,mr)
% ======================================================================= %
% Function to read header and body of a RINEX-file
% ----------------------------------------------------------------------- %
% Input:
%   filename ... RINEX-filename                                 string
%   projectname ... project foldername                          string
%   maxSat ... maximum number of satellites per GNSS            int[1x1]
%   mr...Master or Rover: 1 or 2                                int[1x1]
%
% Output: 
%   header ... header informations                              struct
%   obs ... observations: gps, glonass, galileo [maxSatx15]     struct
%   date ... year, month, day, hour, min, sec                   double[nx6]
%   sow ... seconds of week                                     double[nx1]
%
% ----------------------------------------------------------------------- %
%   @author(s): F. Zimmermann, C. Eling
%   @last modified: 23.03.2017
%   @contact: zimmermann@igg.uni-bonn.de
% ======================================================================= %

% ----- Initialization of output structure and error code --------------- %
errNum = 0;

header = [];
obs = [];
date = [];
sow = [];

% ----- Check function input -------------------------------------------- %
if nargin~=4
    errNum = 1003;
    return
end

if ~ischar(filename) || ~ischar(projectname)
    errNum = 1004;
    return
end

if mr<1 || mr>2
    errNum = 1005;
    return
end


% ----- Read RINEX header ----------------------------------------------- %
[errNum, header] = readRinexHeader(filename);

% ----- Get number of observation epochs -------------------------------- %
% date of first observation
yearF = header.time{1,1};
monthF = header.time{1,2};
dayF = header.time{1,3};
hourF = header.time{1,4};
minF = header.time{1,5};
secF = header.time{1,6};
yearNumF = str2double(yearF);
monthNumF = str2double(monthF);
dayNumF = str2double(dayF);
hourNumF = str2double(hourF);
minNumF = str2double(minF);
secNumF = str2double(secF);

% date of last observation
yearL = header.time{2,1};
monthL = header.time{2,2};
dayL = header.time{2,3};
hourL = header.time{2,4};
minL = header.time{2,5};
secL = header.time{2,6};
yearNumL = str2double(yearL);
monthNumL = str2double(monthL);
dayNumL = str2double(dayL);
hourNumL = str2double(hourL);
minNumL = str2double(minL);
secNumL = str2double(secL);

% get first two observation epochs from RINEX file
interval = zeros(2,6);
count = 1;
fid = fopen(filename, 'r');
dummy = textscan(fid,'%s ',1,'headerlines',header.end-1);
while count<3
    % get next line of file
    tline = fgetl(fid);
    
    yy = str2double(tline(2:3));
    mm = str2double(tline(5:6));
    dd = str2double(tline(8:9));
    
    cmpDate = yy==(yearNumF-2000) & mm==monthNumF & dd==dayNumF;
    
    if cmpDate 
        int1 = textscan(tline,'%n %n %n %n %n %n %*n %*38c',1);
        interval(count,:) = cell2mat(int1);
        interval(count,1) = interval(count,1)+2000;
                
        % increment counter
        count = count + 1;
    end
end
fclose(fid);
% convert observation epoch to SOW
tmp1 = getSOW(interval(1,:)');
tmp2 = getSOW(interval(2,:)');

% compute complete observation duration in SOW
firstEpoch = getSOW([yearNumF,monthNumF,dayNumF,hourNumF,minNumF,secNumF]);
lastEpoch = getSOW([yearNumL,monthNumL,dayNumL,hourNumL,minNumL,secNumL]);
duration = lastEpoch-firstEpoch;

% compute number of observation epochs
numEpochs = duration/(tmp2-tmp1);

% ----- Read RINEX observations ----------------------------------------- %
% initialize output structure
obs = struct('obsG',cell(numEpochs,1),'obsR',cell(numEpochs,1),...
             'obsE',cell(numEpochs,1),'satG',cell(numEpochs,1),...
             'satR',cell(numEpochs,1),'satE',cell(numEpochs,1));
date = zeros(numEpochs,6);
sow = zeros(numEpochs,1);
         
% initialize current number of observation epoch
epochNum = 1;

% move through file and look for observation times
fid = fopen(filename, 'r');
dummy = textscan(fid,'%s ',1,'headerlines',header.end-1);
while 1
    
    if mod(epochNum,250)==0
        disp(['      observation epoch: ' num2str(epochNum)])
    end
    
    % get next line of file
    tline = fgetl(fid);
    % check for end of file
    if ~ischar(tline)
        break
    end
    
    yy = str2double(tline(2:3));
    mm = str2double(tline(5:6));
    dd = str2double(tline(8:9));
    
    cmpDate1 = yy==(yearNumF-2000) & mm==monthNumF & dd==dayNumF;
    cmpDate2 = yy==(yearNumL-2000) & mm==monthNumL & dd==dayNumL;
    if cmpDate1 || cmpDate2
        
        % ----- identify date and satellite PRNs ------------------------ %
        %epoch1 = textscan(tline,'%n %n %n %n %n %n %n %38c',1);
        epoch1 = textscan(tline,'%n %n %n %n %n %n %n %[^c]',1);
        line1 = epoch1{1,8};
        line1 = line1{1};
        % check number of observed satellites
        if line1(1)=='1' || line1(1)=='2' || line1(1)=='3'
            if isnan(str2double(line1(2))) == 1
                nSat = str2double(line1(1));
            else
                nSat = str2double([line1(1),line1(2)]);
            end
        else if line1(1)=='0'
             nSat = str2double(line1(1:2));
            else
             nSat = str2double(line1(1));
            end
        end
        
        % if number of observed satellites >12, PRNs stored in two lines
        % epoch1 -> line 1, epoch2 -> line 2
        if nSat>12
            epoch2 = textscan(fid,'%[^\n]',1);
            line2 = cell2mat(epoch2{1,1});
            dummy = fgets(fid);
        else
            epoch2 = [];
        end
        % if number of observed satellites >24, PRNs stored in three lines
        % epoch1 -> line 1, epoch2 -> line 2, epoch3 -> line 3
        if nSat>24
            epoch3 = textscan(fid,'%[^\n]',1);
            line3 = cell2mat(epoch3{1,1});
            dummy = fgets(fid);
        else
            epoch3 = [];
        end

        % ----- get observation epoch ----------------------------------- %
        for r = 1:6
            date(epochNum,r) = epoch1{1,r};
        end
        date(epochNum,1) = date(epochNum,1)+2000;
        
        % convert observation epoch to SOW
        sow(epochNum,1) = getSOW(date(epochNum,:)');
        
        % ----- get satellite PRNs -------------------------------------- %
        if ~isempty(epoch2) && ~isempty(epoch3)
            prnStr = [line1,line2,line3];
        elseif ~isempty(epoch2)
            prnStr = [line1,line2];
        else
            prnStr = line1;
        end
        tf = isletter(prnStr);
        satTypes = prnStr(tf);
        satTypeID = zeros(nSat,1);
        satPrn = zeros(nSat,1);
        prnIdx = find(tf);
        for k = 1:length(satTypes)-1
            satPrn(k) = str2double(prnStr(prnIdx(k)+1:prnIdx(k+1)-1));
            if strcmp(satTypes(k),'G')
                satTypeID(k) = 1;
            elseif strcmp(satTypes(k),'R')
                satTypeID(k) = 2;
            elseif strcmp(satTypes(k),'E')
                satTypeID(k) = 3;
            end
        end
        satPrn(end) = str2double(prnStr(prnIdx(end)+1:end));
        if strcmp(satTypes(end),'G')
            satTypeID(end) = 1;
        elseif strcmp(satTypes(end),'R')
            satTypeID(end) = 2;
        elseif strcmp(satTypes(end),'E')
            satTypeID(end) = 3;
        end
        
        obs(epochNum,1).satG = satPrn(satTypeID == 1);
        obs(epochNum,1).satR = satPrn(satTypeID == 2);
        obs(epochNum,1).satE = satPrn(satTypeID == 3);
               
        % ----- get observations ---------------------------------------- %
        [errNum,obsOut] = readRinexObs(fid,header,nSat);
        
        % ----- allocate observations to output file -------------------- %
        gps = zeros(maxSat,15);
        glo = zeros(maxSat,15);
        gal = zeros(maxSat,15);

        ind = find(header.obsCol(header.obsCol~=0),1,'last'); 
        if header.snrMapped == 1
            ind = ind+1;
            gps(obs(epochNum,1).satG,1:ind) = obsOut(satTypeID==1,:);
            glo(obs(epochNum,1).satR,1:ind) = obsOut(satTypeID==2,:);
            gal(obs(epochNum,1).satE,1:ind) = obsOut(satTypeID==3,:);
        else
            % ind = ind+1
            % gps(obs(epochNum,1).satG,1:ind) = obsOut(satTypeID==1,:);
            % glo(obs(epochNum,1).satR,1:ind) = obsOut(satTypeID==2,:);
            % gal(obs(epochNum,1).satE,1:ind) = obsOut(satTypeID==3,:);
            gps(obs(epochNum,1).satG,1:ind) = obsOut(satTypeID==1,header.obsCol(header.obsCol~=0)');
            glo(obs(epochNum,1).satR,1:ind) = obsOut(satTypeID==2,header.obsCol(header.obsCol~=0)');
            gal(obs(epochNum,1).satE,1:ind) = obsOut(satTypeID==3,header.obsCol(header.obsCol~=0)');
        end
        
        obs(epochNum,1).obsG = gps;
        obs(epochNum,1).obsR = glo;
        obs(epochNum,1).obsE = gal;
              
        % increment current number of observation epoch
        epochNum = epochNum+1;
    end
   
end
fclose(fid);

% find and delete epochs with missing data
missEpochs = find(sow==0);
sow(missEpochs) = [];
date(missEpochs,:) = [];
obs(missEpochs) = [];

return;