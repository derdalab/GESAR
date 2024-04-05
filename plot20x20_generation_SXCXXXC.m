%close all;

addpath (genpath(pwd));

D = '';
F = 'VT_unfiltered_Feb.txt';
SAVEfile = 'C3C.txt';  % keep blank if don't want to save


SET{1} = 1:3;         % B-SX4, captured: size 6e7 PFU
SET{2} = 7:9;         % SX4, captured: size 2e4 PFU
SET{3} = 4:6;         % B-SX4, captured on biotin-blocked beads: 3e3 PFU
SET{4} = 10:12;         % input of B-SX4,
SET{5} = 16:18;         % input of SX74,
SET{6} = 13:15;         % input of B-SX4, captured on biotin-blocked beads

PFU(1) = 6e7;
PFU(2) = 2e4;
PFU(3) = 3e3;
PFU(4) = 1e10;
PFU(5) = 1e10;
PFU(6) = 1e10;

INPUTS = [4 5 6];

TEST_SET = 1;

SHOWaminoACIDS = [19 20 21 22 23 24];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SAVEDir = '';
saveIMG = 'YEBreactiveC3C';
SAVEto = 'YEBreactiveC3C.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% find the largest element in the SET %%%%%%%%%%%%%%%%%%%
MAX = 0;
for i=1:numel(SET)
    if max(SET{i}) > MAX
        MAX = max(SET{i});
    end
end

% define the numeric values for the input columns
input_columns = [];
for i = 1:numel(INPUTS)
    input_columns = [input_columns SET{INPUTS(i)}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REREAD = 1; % change to zero if you're rerun the script 

if REREAD
    [Nu0, AA0, Fr0] = readMulticolumn('Dir', D, 'File', F, ...
                                    'column', 1:MAX,...
                                    'skip', 2, 'output', 'raw');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% find the most common SDBs
cNuc = char(Nu0);
SDB = cellstr( cNuc(:,7:27) );

[uSDB,frSDB]=uniqueCOMB( SDB, ones (size(SDB,1), 1 ) );

% the barcodes of C3C are 2nd, 3rd, 4th and 5th most common; add them up
SDB_numbers= 2:5;

C3C    = 'TCT~~~TGT~~~~~~~~~TGTGGTGGAGGT'    ;

% found all the library members that have this SDB
S1= [];
for j = SDB_numbers
    
    SDB{j} = ['~~~~~~' uSDB{j} '~~~~~~~~~~~~~~~~~~~~~~~~'];
    [temp, ~] = unfiltered2filtered(Nu0, SDB{j}, 1); fprintf('.');
    S1 = [S1; temp];
end

% trim the library to that SDB
trNuc = cNuc(S1, numel(SDB{j})+1 : end);


% filter out the SX4 library
[S2, ~] = unfiltered2filtered(trNuc,C3C, 1); 

Nu = Nu0 ( S1(S2), :);
AA = AA0 ( S1(S2), :);
Fr = Fr0 ( S1(S2), :);
%%
% select only the aminoacids you want to see, SELECT ONLY SX4 motifs!!
cAA = char(AA);
AA=cellstr(cAA(:,SHOWaminoACIDS));

%%%%%% combine the frequencies

[AA,Fr]=uniqueCOMB(AA,Fr);

%%%% select only the amino acids reliably present in the input
inFr = sum( Fr(: ,input_columns), 2);
reliable = find( inFr > 2);

%% plot the input
               
VAR =       {'Nuc',AA,...
                   'LookAt',[1 3 4 5],...
                   'disp400x400',1,'disp20x400',0,...
                   'maxFreq', 10000, 'scale', 'log2',...
                   'grid',1,...
                   'hydrophobicity','janin',...
                   'sample',0};

visual400x400gen('figure', figure(1), 'Fr', inFr, VAR{:});
visual400x400gen('figure', figure(2), 'Fr', sum( Fr(:, SET{1}), 2), VAR{:});
visual400x400gen('figure', figure(3), 'Fr', sum( Fr(:, SET{2}), 2), VAR{:});
visual400x400gen('figure', figure(4), 'Fr', sum( Fr(:, SET{3}), 2), VAR{:});

Nu = Nu ( reliable, :);
AA = AA ( reliable, :);
Fr1 = Fr ( reliable, :);
inFr = inFr(reliable, :);
 %%                    
for i =1:6
    Sa(i) = PFU(i) / sum(sum(Fr(:,SET{i})));
    Si{i} = Sa(i)*sum(Fr1(:,SET{i}),2);
end

VAR =       {'Nuc',AA,...
                   'LookAt',[1 3 4 5],...
                   'disp400x400',1,'disp20x400',0,...
                   'maxFreq', 3000, 'scale', 'log2',...
                   'grid',1,...
                   'hydrophobicity','janin',...
                   'sample',0};

visual400x400gen('figure', figure(5), 'Fr', inFr, VAR{:});
visual400x400gen('figure', figure(6), 'Fr', sum( Fr1(:, SET{1}), 2), VAR{:});
visual400x400gen('figure', figure(7), 'Fr', sum( Fr1(:, SET{2}), 2), VAR{:});
visual400x400gen('figure', figure(8), 'Fr', sum( Fr1(:, SET{3}), 2), VAR{:});

                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%
C = (Si{1} - Si{2} - Si{3}) ./ ( Si{4} + Si{5} + Si{6} );
%C2 = (Si{1} ) ./ ( Si{4} + Si{5} + Si{6} );

%find inities in the ratio and replaced them by a finite large number.
INF = find(isinf(C));
MAX = max(C(C<inf));
C(INF,:)=MAX*2;

%zero doesn't plot well on log-scale; replace zeros by some small numbers
NEG = find(C(:,1)<0);
ZER = find(C(:,1)==0);
%ZER2 = find(C2(:,1)==0);

% define small number as the smallest non-zero number in the set
%make them the same as in SX4 set
%MIN=1.0895e-04;
MIN=2^-13;

% MIN = min(C(C>0));
% MIN2 = min(C2(C2>0));
%%%%%%% FIX THIS LATER TO MAKE MORE INTELLIGENT CONNECTION BETWEEN SETS

% replace zeros & negative numbers by a 1/2 and 1/3 of the smallest number
C(ZER,:)=MIN/3;
C(NEG,:)=MIN/4;

% find zero/zero = NaN and replace them by 1/5 of the smallest number
NANS = find(isnan(C(:,1)));
C(NANS,:)=MIN/5;

% note that in the "reliable dataset, there will be no INF and no NANS

%%


% data cannot be sampled because the inputs are floating point ratios not integer counts
[dist400]=visual400x400gen('Nuc',AA, 'Fr',log2(C)-log2(MIN/4),...
                           'LookAt',[1 3 4 5],...
                           'disp400x400',1,'disp20x400',0,...
                           'maxFreq', 10,... log2(MAX*2)-log2(MIN/4),...
                           'scale', 'lin',...
                           'grid',1,...
                           'figure', figure(10),...
                           'hydrophobicity','janin',...
                           'sample',0);
                       saveas(gcf,saveIMG,'epsc')


S = char( 32*ones(size(AA,1),1) );
toSave = [char(AA)     S num2str(C,'%10.5e')  ];

if ~isempty(SAVEfile)
    fs = fopen(fullfile(SAVEDir,SAVEfile),'w');
    RET = char(10*ones(size(toSave,1),1));
    fprintf( fs, '%s\r\n', [toSave RET]');
    fclose all;
end


%%%%%% plot histogram of ratio
% figure(30);
% 
% [b,x]=hist( log2(C)-log2(MIN/4), 0:0.5: (log2(MAX*2)-log2(MIN/4)) );
% bar(x,b);
% set(gca,'yscale','log',...
%        'xtick', 0 : 1 : ceil(log2(MAX*2)-log2(MIN/4)) ,...
%        'TickDir','out');
% xlim([0; ceil(log2(MAX*2)-log2(MIN/4)) ]);
% saveas(gcf,[saveIMG 'Hist'],'epsc')
% 
% % figure(40);
% % 
% % [b,x]=hist( log2(C2)-log2(MIN2/4), 0:0.5: (log2(MAX*2)-log2(MIN2/4)) );
% % bar(x,b);
% % set(gca,'yscale','log',...
% %        'xtick', 0 : 1 : ceil(log2(MAX*2)-log2(MIN2/4)) ,...
% %        'TickDir','out');
% % xlim([0; ceil(log2(MAX*2)-log2(MIN2/4)) ]);
% % saveas(gcf,[saveIMG 'Hist'],'epsc')

figure(30);


[b,x]=hist( log2(C), floor(log2(MIN/4)) :0.5: log2(MAX*2) );
bar(x,b);
  
   set(gca,'yscale','log',...
       'xtick', floor(log2(MIN/4)) : 1 : ceil(log2(MAX*2)) ,...
       'TickDir','out');

%xlim([log2(MIN/5); ceil(log2(MAX*2)) ]);
xlim([log2(MIN/5); -4 ]);
saveas(gcf,[saveIMG 'Hist'],'epsc')

figure(100); 
plot(Si{1},Si{2},'.k')
hold on
plot(Si{1},Si{3},'.r')
plot(Si{1},Si{4},'.c')
LIM = [0.001 1e8];
plot(LIM, LIM,'-b');
set(gca,'xscale','log','yscale','log');
xlim(LIM);
ylim(LIM);


