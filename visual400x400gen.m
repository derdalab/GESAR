function [out400,out20,Nuc4, Fr4] = visual400x400gen(varargin)
 
% by default display 400x400 plot but not 20x400
disp400x400= 1;
disp20x400 = 0;
scale='lin';

BORDER = 0.7*[1 1 1];
ZERO   = 0.1*[1 1 1];

% by default minimum frequency is 0
minFreq = 0;

% if you want to sample the library, change this 
SAMPLE=0;
 
% default amino acids you want to look at
LOOKAT = 1:4;

% % normalize by codon frequency or convert to PPM
NORM  = 0;
seq   = 'aa';

% order in which the amino acids appear
hydrophobicity = 'alphabetical'; 

% default colormap
COLORMAP = jet;

% default figure handle
h = 10;

if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'nuc',         Nuc = varargin{ni+1};
            case 'fr',          Fr = varargin{ni+1};
            case 'seq',         seq=varargin{ni+1};
            case 'lookat',      LOOKAT = varargin{ni+1};
            case 'disp400x400', disp400x400 = varargin{ni+1}; 
            case 'disp20x400',  disp20x400 = varargin{ni+1}; 
            case 'sample',      SAMPLE = varargin{ni+1}; 
            case 'scale',       scale = varargin{ni+1};
            case 'maxfreq',     maxFreq = varargin{ni+1};
            case 'minfreq',     minFreq = varargin{ni+1};
            case 'norm',        NORM = varargin{ni+1};
            case 'grid',        GRID = varargin{ni+1};
            case 'hydrophobicity', hydrophobicity = varargin{ni+1};
            case 'colormap',    COLORMAP = varargin{ni+1};
            case 'figure',      h = varargin{ni+1};
        end
    end
end 

%%%%%%%%%%% the order in which amion acids appear on 400x400 plot %%%%%%%
if strcmp(hydrophobicity,'alphabetical')    
    AA = [  'A';'C';'D';'E';'F';'G';'H';'I';'K';
            'L';'M';'N';'P';'Q';'R';'S';'T';'V';'W';'Y'];
elseif strcmp(hydrophobicity,'janin')   
    AA = [  'R';'K';'Q';'E';'D';'N';'Y';'P';'T';
            'S';'H';'A';'G';'W';'M';'F';'L';'V';'I';'C'];
else
    AA = [  'A';'C';'D';'E';'F';'G';'H';'I';'K';
            'L';'M';'N';'P';'Q';'R';'S';'T';'V';'W';'Y'];
end

if SAMPLE
    % take a random subset of N=SAMPLE sequences. The sum has to be=SAMPLE 
    Fr = SA(Fr,SAMPLE);
    disp(['sampled ' num2str(sum(Fr)) ' total reads']);
end

NN = numel(AA);

% combine non-unique reads; recalculate number of unique reads
[Nuc,Fr] = uniqueCOMB(Nuc,Fr);

disp(['loaded ' num2str(numel(Fr)) ' unique and ' ...
      num2str(sum(Fr)) ' total reads']);

% convert the sequences to amino acids if nucleotides are suplied
if strcmp(seq,'nt')
    cNuc = nt2aacell(Nuc,1);
else
    cNuc = char(Nuc);
end

% calculate the truncated frequencies
disp('calculating the 4-mer frequencies');

% extract relevant positions
cNuc4 = cNuc(:,LOOKAT);
Nuc4 = cellstr(cNuc4);

[Nuc4,Fr4] = uniqueCOMB(Nuc4,Fr);
cNuc4 = char(Nuc4);
disp('done');


% calculate the frequencies in the 400x400 matrix
disp('calculating the frequency in 400x400 sequence matrix');

temp2 = cell2numaa(cNuc4, Fr4, AA);

%%%%%%%%%%%%% end of calculate the frequencie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NORM
% normalize the frequencies using parts per million
    temp2 = 10^6 * temp2./ sum(Fr);
end
disp('done');



%%%%% add ligh grey grid to 400x400 plot
% values added are -1 (easy to find later)
if GRID
    
    temp3=[];
    
    Xline = -ones(NN^2,1);
    start=1;
    
    for i=1:NN
        temp3 = [temp3  Xline temp2(: , start:start+NN-1)  ];
        start=start+NN;
    end
    
    Yline = -ones(1,NN^2+NN);
    start=1;
    temp4=[];
    
    for i=1:NN
        temp4 = [temp4;  Yline; temp3(start:start+NN-1, : )  ];
        start=start+NN;
    end
    
    out400=temp4;
else
    out400=temp2;
end

% find zeros 
ZEROS = find(out400==0);
disp(['There are ' num2str(numel(ZEROS)) 'zeros' ]);
LINES = find (out400 ==  -1);

if disp400x400   
    figure(h);
    disp('plotting the data');
        
    if strcmp(scale,'log')

        temp3 = log10(out400);
        % replace all values in negative and zero elements
        if GRID
            temp3(LINES) = -log10(maxFreq)/31;
            temp3(ZEROS) = -log10(maxFreq)/62;
        end
        
        plotting=temp3;

        imagesc(plotting,[-log10(maxFreq)/31 log10(maxFreq) ]);
        % color range
    elseif strcmp(scale,'log2')

        temp3 = log2(out400);
        % replace all values in negative and zero elements
        if GRID
            temp3(LINES) = -log2(maxFreq)/31;
            temp3(ZEROS) = -log2(maxFreq)/63;
        end
        
        plotting=temp3;
        imagesc(plotting,[-log2(maxFreq)/31 log2(maxFreq) ]);
        % color range       
    else 
        temp3=out400;
        if GRID
            temp3(LINES) = -(1/31)*maxFreq;
            temp3(ZEROS) = -(1/63)*maxFreq;
        end
        plotting = temp3; 

        imagesc(plotting,[-(1/31)*maxFreq maxFreq ]);

    end

    jet2=COLORMAP;

    % this part is pretty sensitive ot overall range of copy numbers
    % you have to be carefull changin this part becaue you can turn low copy 
    % numbers to "black" by accident.
    
    jet2(1,:)= BORDER;  % grey color for a border
    jet2(2,:)= ZERO;

    colormap(jet2);
    colorbar;

    if GRID
        set(gca,'YTick', (1:(NN+1):(NN+1)^2-1)+NN/2, 'YTickLabel',AA,'TickDir','out');
        set(gca,'XTick', (1:(NN+1):(NN+1)^2-1)+NN/2, 'XTickLabel',AA);
    else
        set(gca,'YTick', 1:NN:NN^2-1, 'YTickLabel',AA,'TickDir','out');
        set(gca,'XTick', 1:NN:NN^2-1, 'XTickLabel',AA);
    end
    drawnow;
end


% this is where 400x400 plot is condensed to 20x400 plot
A20x400=[];
 
for jj=1:20
    A20x400=[A20x400; sum(temp2(jj:20:400,:),1)];
end

out20=A20x400;

if disp20x400
    figure(11);
    plotting2 = (A20x400);
    imagesc(plotting2,[0 50]);

    jet2=COLORMAP;
    jet2(1,:)= 0.1*[1 1 1];

    colormap(jet2);
    colorbar;
    set(gca,'YTick', 1:1:20, 'YTickLabel',AA,'TickDir','out');
    set(gca,'XTick', 1:20:400-1, 'XTickLabel',AA);
    drawnow;
end


end
