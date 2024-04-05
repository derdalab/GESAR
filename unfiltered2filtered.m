function [OUT,ISCODON] = unfiltered2filtered(IN, filter, h, varargin)
% this function take in a cell array of strings that contains nucleotide
% sequences (IN variable) and also a filter, in form of 
% filter = 'TCT~~~TGT~~~~~~~~~TGTGGTGGAGGT';
% where ~ designates any nucleotide and the other nucleotides are part of
% filter. It then look for all nucleotides that are withing a h -distance
% of the filter (use h=1 for loose filtering and h=0 for super sctrict
% filtering)
% usage for filtering of true-blue CxCxxxCGGG reads from Nu cell array
% [OUT] = unfiltered2filtered(Nu, 'TCT~~~TGT~~~~~~~~~TGTGGTGGAGGT', 0);  
% the output OUT is an index of the reads the pass the filter
% ISCODON, is a yes no answer whether codon from the specific position
% comes from the allowed set. If the set is not specificied, the program
% has internally stored TriNuc set by default
% TriNuc = ['(TGG)|(GGT)|(CTG)|(CCG)|(CGT)|(AAC)|(AAA)|(ATG)|(GAA)|(CAG)' ...
%          '|(GCT)|(CAT)|(TCT)|(GAC)|(ACT)|(TTC)|(GTT)|(ATC)|(TAC)'];

codons = {'TGG', 'GGT', 'CTG', 'CCG', 'CGT', 'AAC', 'AAA', 'ATG',...
          'GAA', 'CAG', 'GCT', 'CAT', 'TCT', 'GAC', 'ACT', 'TTC',...
          'GTT', 'ATC', 'TAC'};

positions = {4:6, 10:12, 13:15, 16:18};
START = 1;
look4codons = 1;

if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'codons',      codons = varargin{ni+1};
            case 'positions',   positions = varargin{ni+1};
            case 'start',       START = varargin{ni+1};
            case 'look4codons', look4codons = varargin{ni+1};
        end
    end
end 

cIN = char(IN);

MAX = min (abs('ATCG '-'~~~~~'));  % calculate maximum allowable difference
Afilter = makearray(filter,size(IN,1));  %turn the filter into array form

OUT = [];

if isnumeric(START)
    
    % truncate the reads to filter length, by default start from position 1
    cIN0 = cIN(:, START:numel(filter)+START-1); 

    
    D = abs(cIN0-Afilter);  % calculate the absolute difference

    % if all values are zero, the sum iz zero, if there is one mismatch with
    % the filter, the sum is 1, if there is more than one, then the read is not
    % considered

    sumD = sum ( ~( (D >= MAX) | (D == 0) ),2 );  

    % the function returns the indices of the reads that match the filter
    OUT = find(sumD<=h);
else
    for i=1: (numel(cIN(1,:)) - numel(filter))+1
        START = i;
        cIN0 = cIN(:, START:numel(filter)+START-1); 

        D = abs(cIN0-Afilter);  % calculate the absolute difference

        sumD = sum ( ~( (D >= MAX) | (D == 0) ),2 );  

    % the function returns the indices of the reads that match the filter
        OUT = union(OUT,find(sumD<=h));
    end
end

if ~look4codons
    ISCODON = [];
    return 
end
%%%%%%%%%%%%%%%%%%%%% this part matches the codons %%%%%%%%%%%%%%%%%%%%%%%
% replicated the codons array to the size of the OUT variable
Cfilter={};
for i=1:numel(codons)
    Cfilter{i} = makearray(codons{i},size(OUT,1));
end

%go through each position 
temp = zeros( size(OUT,1), numel(positions));

for i=1:numel(positions)
    for j=1:numel(Cfilter)
        temp(:,i) = temp(:,i) + ~sum( cIN(OUT,positions{i})~=Cfilter{j}, 2);
    end
end

ISCODON = temp;

end

function [aSEQ] = makearray(SEQ,SIZE)

    aSEQ = ones(SIZE,numel(SEQ));
    
    for i=1:numel(SEQ)
        aSEQ(:,i)=SEQ(i);
    end
    
    aSEQ = char(aSEQ);

end