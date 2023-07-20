function [errNum, header] = iggReadRinexHeader(filename)
% ======================================================================= %
% Function to read the header of a RINEX-file
%
%   based on: 
%       'RINEX: The Receiver Independent Exchange Format Version 2.11'
% ----------------------------------------------------------------------- %
% Input:
%   filename ... RINEX filename                                 str[-]
%
% Output: 
%   errNum ... error code                                       int[1x1]
%   header ... RINEX header informations                        struct[1x1]
%   field description:
%   ->    marker ... marker name                                str[-]
%   ->    recNum ... receiver serial number                     str[-]
%   ->   recType ... receiver type                              str[-]
%   ->    antNum ... antenna serial number                      str[-]
%   ->   antType ... antenna type                               str[-]
%   ->   posCode ... approximate receiver position [m] (X,Y,Z)  double[1x3]
%   ->  deltaAnt ... antenna offset [m] (X,Y,Z)                 double[1x3]
%   -> nObsTypes ... number of observation types                int[1x1]
%   ->  obsTypes ... abbreviations of observation types         str[-]
%   ->    obsCol ... columns of observation types               int[15x1]
%         -> C1,L1,S1,P2,L2,S2,C5,L5,S5,C7,L7,S7,C8,L8,S8       int[1x1]
%   ->    obsCol ... columns of observation types               struct[1x1]
%   -> snrMapped ... indicator for mapped SNR values            int[1x1]
%   ->      time ... first/last obs times [str] (y,m,d,h,m,s)   cell[2x6]
%   ->   leapSec ... number of leap seconds                     int[1x1]
%   ->       end ... number of last header line                 int[1x1]
% ----------------------------------------------------------------------- %
%   @author(s): Florian Zimmermann
%   @last modified: 11.11.2016
%   @contact: zimmermann@igg.uni-bonn.de
% ======================================================================= %

% ----- Initialization of output structure and error code --------------- %
% init error return
errNum = 0;
% output structure
header = struct('marker',[],'recNum',[],'recType',[],'antNum',[],...
                'antType',[],'posCode',[],'deltaAnt',[],'nObsTypes',[],...
                'obsTypes',[],'obsCol',[],'snrMapped',[],'time',[],...
                'leapSec',[],'end',[]);
header.time = cell(2,6);

% initialize counter
count = 1;

% ----- Check function input -------------------------------------------- %
if nargin~=1
    errNum = 1001;
    return
end

if ~ischar(filename)
    errNum = 1002;
    return
end

% ----- Read RINEX header ----------------------------------------------- %
% open file for reading
fid = fopen(filename, 'r');
while 1
    % get next line of file
    tline = fgetl(fid);
    
    % find marker name
    if (strfind(tline,'MARKER NAME'))
        marker = textscan(tline,'%s',1);
        header.marker = marker{1,1}{1,1};
    end
    
    % find number and type of receiver
    if (strfind(tline,'REC #'))
        recNumType = textscan(tline,'%20c %20c',1);
        header.recNum = recNumType{1,1};
        header.recType = recNumType{1,2};
    end
    
    % find number and type of antenna
    if (strfind(tline,'ANT #'))
        antNumType = textscan(tline,'%20c %20c',1);
        header.antNum = antNumType{1,1};
        header.antType = antNumType{1,2};
    end
        
    % find approximate receiver position
    if (strfind(tline,'APPROX POSITION XYZ'))
        pos_code = textscan(tline,'%n %n %n ',1);
        header.posCode = cell2mat(pos_code);
    end
    
    % find antenna offset
    if (strfind(tline,'ANTENNA: DELTA H/E/N'))
        delta_ant = textscan(tline,'%n %n %n ',1);
        header.deltaAnt = cell2mat(delta_ant);
    end
    
    % find observation types
    if (strfind(tline,'# / TYPES OF OBSERV'))
        obsTypeLine1 = textscan(tline,'%n %[^#]');
        header.nObsTypes = obsTypeLine1{1,1};
        obsTypes1 = obsTypeLine1{1,2}{1,1};    
        
        % check for second line of observation types
        if header.nObsTypes>9
            count = count + 1;
            % get next line of file
            tline = fgetl(fid);
            obsTypeLine2 = textscan(tline,'%[^#]');
            obsTypes2 = obsTypeLine2{1,1}{1,1};
        else
            obsTypes2 = [];
        end
        
        % merge observation types
        obsTypes = [obsTypes1 ' ' obsTypes2];
        header.obsTypes = obsTypes;
        
        % get sorting of observation types in RINEX file
        obsCol = zeros(15,1);            
        sf = ~isspace(header.obsTypes);
        types = header.obsTypes(sf);
        for m = 1:header.nObsTypes
            tmp = types(m*2-1:m*2);
            if strcmp(tmp,'C1')
                obsCol(1) = m;
            elseif strcmp(types(m*2-1:m*2),'L1')
                obsCol(2) = m;
            elseif strcmp(types(m*2-1:m*2),'S1')
                obsCol(3) = m;
            elseif strcmp(types(m*2-1:m*2),'P2')
                obsCol(4) = m;
            elseif strcmp(types(m*2-1:m*2),'L2')
                obsCol(5) = m;
            elseif strcmp(types(m*2-1:m*2),'S2')
                obsCol(6) = m;
            elseif strcmp(types(m*2-1:m*2),'C5')
                obsCol(7) = m;
            elseif strcmp(types(m*2-1:m*2),'L5')
                obsCol(8) = m;
            elseif strcmp(types(m*2-1:m*2),'S5')
                obsCol(9) = m;    
            elseif strcmp(types(m*2-1:m*2),'C7')
                obsCol(10) = m;
            elseif strcmp(types(m*2-1:m*2),'L7')
                obsCol(11) = m;
            elseif strcmp(types(m*2-1:m*2),'S7')
                obsCol(12) = m;
            elseif strcmp(types(m*2-1:m*2),'C8')
                obsCol(13) = m;
            elseif strcmp(types(m*2-1:m*2),'L8')
                obsCol(14) = m;
            elseif strcmp(types(m*2-1:m*2),'S8')
                obsCol(15) = m;
            end
        end
        header.obsCol = obsCol;
        
        % check for mapped SNR values
        if header.obsCol(3)==0
            header.snrMapped = 1;
        else
            header.snrMapped = 0;
        end
    end
        
    % find time of first observation
    if (strfind(tline,'TIME OF FIRST OBS'))
        start_obs = textscan(tline,'%s %s %s %s %s %s %*s %*s %*s %*s  ',1);
        year = start_obs{1};
        month = start_obs{2};
        day = start_obs{3};
        hour = start_obs{4};
        min = start_obs{5};
        sec = start_obs{6};
        header.time(1,:) = [year,month,day,hour,min,sec]; 
    end
    
    % find time of last observation
    if (strfind(tline,'TIME OF LAST OBS'))
        last_obs = textscan(tline,'%s %s %s %s %s %s %*s %*s %*s %*s  ',1);
        year = last_obs{1};
        month = last_obs{2};
        day = last_obs{3};
        hour = last_obs{4};
        min = last_obs{5};
        sec = last_obs{6};
        header.time(2,:) = [year,month,day,hour,min,sec]; 
    end
    
    % find current number of leap seconds
    if (strfind(tline,'LEAP SECONDS'))
        leapS = textscan(tline,'%n',1);
        header.leapSec = cell2mat(leapS);
    end
    
    % find end of header
    if (strfind(tline,'END OF HEADER'))
        header.end = count;
        break;
    end
    
    % increment counter
    count = count + 1;
end
% close file
fclose(fid);

return;