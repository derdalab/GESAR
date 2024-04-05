function [uSeqLo, sumNum, Occur] = uniqueCOMB (SeqLo, NumLo, varargin)
% it's just like an ordinary uniqueCOMB function that takes SeqLo | NumLo
% AAATA 3
% AAATA 2
% AAAAA 4
% and converts to 
% AAATA 5
% AAAAA 4
% but here SeqLo is an array of NxM instead of Nx1
% The script also has new input; usage: 
% [uNu, sFr] = uniqueCOMB(Nu, Fr, 'sortby', [2 3 4])

sortby = 1:size( NumLo, 2); 
% For Nx1 case, sortby = 1; for NxM case, it will sort by sum of columns

    if exist('varargin','var')
        L = length(varargin);
        if rem(L,2) ~= 0 
            error('Parameters/Values must come in pairs.'); 
        end

        % read input variables
        for ni = 1:2:L
            switch lower(varargin{ni})
                case 'sortby',      sortby = varargin{ni+1};
                case 'other var',   dummy_var = varargin{ni+1};   
            end
        end
    end

% if input is not cell, convert it to cell
if ~iscell(SeqLo)
    cellstr(SeqLo);
end

[SeqLo,IX] = sort(SeqLo);  
NumLo = NumLo(IX,:);   % sort both arrays; 
% the line above used to be NumLo = NumLo(IX) for 1xN array.
clear IX

[uSeqLo, M, ~] = unique (SeqLo, 'first');  
M2 = [M(2:end); numel(SeqLo)+1] ;          % M2 is the last element +1
Occur = M2-M; 
clear M2;

% for truncated seqeunces that occured only once, 
% the frequencies are the same as those of the original seqeunces

I1 = find(Occur==1); % indices of sequences that occured once
% M(Il) are the indices of the sequences that occured only once in the
% original array SeqLo
% I1 are the indices of the sequences that occured once in the unique array

sumNum = nan(size(uSeqLo,1), size(NumLo,2) );
% make an empty array os NaNs with N2xM size, where N2 is height as unique
% sequence array and M is the original width of the NumLo array)
% Previous version used to say sumNum = nan(size(uSeqLo));

sumNum( I1, : ) = NumLo( M(I1), : );
% previous version used to say sumNum( I1 ) = NumLo( M(I1) );
clear I1;

% if truncations were found more than once, the frequencies must be added

I2 = find(Occur>1);  
Nu = Occur(I2);   % this is the number of times each sub-seqeunce was found

for iii = 1:numel(I2)
    temp = M(I2(iii));
    sumNum( I2(iii) , : ) = sum (NumLo( temp:temp+Nu(iii)-1, : ), 1 );

end

% previous code in the line above 
% sumNum( I2(iii) ) = sum (NumLo( temp:temp+Nu(iii)-1 ) );

% the next three lines is where there is only drematic different between
% Nx1 code and NxM code. The sorting in Nx1 is obvious but in NxM is not;
% either it has to be defined via 'sortby' function or somehow the code has
% to determine which rows were sorted at first and use the same. I think it
% is safer to supply sortby variable as varargin. In Nx1 case, sortby=1
% [sumNum, IX] = sort(sumNum, 'descend');
%  uSeqLo = uSeqLo(IX);
%  Occur  =  Occur(IX);

[~, IX] = sort( sum(sumNum(:,sortby), 2), 'descend');
sumNum =  sumNum(IX,:);
uSeqLo =  uSeqLo(IX);
Occur  =  Occur(IX);

end   