function [Nuc, AA, FrOUT] = readMulticolumn(varargin)
% this function reads a multi column file Nuc AA Fr1 Fr2 Fr3 etc
% usage  
% readMulticolumn('file', filename, 'dir', directoryname,...
%                   'column', [2 4 9],...
%                   'skip', 2,...
%                   'output',output_type)

    output = 'raw';     % other output_type: 'normalized' 'normalized+1'
    skip = 2;           % by default skip the first two rows
    filtered=0;         % by default read unfiltered files with Nuc and AA

    if exist('varargin','var')
        L = length(varargin);
        if rem(L,2) ~= 0 
            error('Parameters/Values must come in pairs.'); 
        end

        % read input variables
        for ni = 1:2:L
            switch lower(varargin{ni})
                case 'file',          	File = varargin{ni+1};
                case 'dir',             Dir=varargin{ni+1};
                case 'column',          column=varargin{ni+1}; 
                case 'skip',            skip=varargin{ni+1}; 
                case 'output',          output=varargin{ni+1};
                case 'filtered',        filtered=varargin{ni+1};
            end
        end
    end
    
    if filtered
        FORM = '%s';
    else
        FORM = '%s %s';
    end

    if strcmp(column,'all')
        % open the file and count the number of columns
        FID = fopen( fullfile(Dir,File),'r');
        for i=1:skip
            temp1 = fgetl(FID); % line(s) to be discarded
        end
        temp1 = fgetl(FID); % useful line
        temp2 = strsplit(temp1,' ');  % split the string
        if isempty(temp2{end})
            column = 1:numel( temp2(3:end-1) );
        else
            column = 1:numel( temp2(3:end) );
        end
        
        fclose(FID);
    end
        
    for i=1:max(column)
        if sum(i==column)
            FORM = [FORM ' %f'];
        else
            FORM = [FORM ' %*f'];
        end
    end

    FORM = [FORM ' %*[^\n]'];

    fh = fopen(fullfile(Dir, File),'r');
    for i=1:skip
        temp=fgetl(fh);
    end

    AllVar = textscan(fh,FORM);
    if filtered
        Nuc = AllVar{1};
        AA  = AllVar{1};
        Fr = cell2mat(AllVar(2:end));
    else        
        Nuc = AllVar{1};
        AA  = AllVar{2};
        Fr = cell2mat(AllVar(3:end));
    end

    clear AllVar;
    fclose all;

    FrOUT = zeros(size(Fr));
    Smax = 0;

    if strcmp(output,'normalized+1')
        % normalize the frequencies and add singleton to all zeros


        for i=1:size(Fr,2)    
            S = sum(Fr(:,i));    
            FrOUT(:,i) = (Fr(:,i) + 1 ) / S;
            Smax = max(S,Smax);
        end
    elseif strcmp(output,'normalized')
        % normalize the frequencies 

        for i=1:size(Fr,2)    
            S = sum(Fr(:,i));    
            FrOUT(:,i) = Fr(:,i) / S;
            Smax = max(S,Smax);
        end
    elseif strcmp(output,'raw')

        FrOUT = Fr;
    else

        error('Output types are "raw", "normalized" or "normalized+1"');
    end
    

end