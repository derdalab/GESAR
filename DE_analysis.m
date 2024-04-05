clear all; close all;
clc
Dir='';
File = 'VT_unfiltered_Feb.txt';
%SAVEto = [File(1:end-4) '']; % keep blank if don't want to save
SET{1} = [19 20 21]; %
SET{2} = [22 23 24]; %
SET{3} = [25 26 27]; %
SET{4} = [28 29 30];
SET{5} = [31 32 33];
SET{6} = [34 35 36];
TEST_SET = 4;
CONTROL_SETS = [1];
HITS2DISPLAY = 20; % maximum numer of hits to display
SHOWaminoACIDS = [1 2 3 4 5];
CLUSTERbyH = 1; % if you want your hits to be clustered by Hamming dist.
PLOT_VOLCANO = 1; % set to 1 if you want to see the actual volcano plot
%%%%%%% volcano plot parameters here
p_cutoff = 0.05; % p-value cutoff (use 0.9 if dont care abt p)
R_cutoff = 5; % ratio cutoff
MaxX=14; % maximum on the X-scale (if plotting volcano)
vert_cutoff = 0.00001; % maximum on the Y-scale (if plotting volcano)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% do not change things beyond this point
%%%%%%%%%%%%%%%%%%%%%% unless you know what you are doing
REREAD=1;
if REREAD
[Nuc, AA, Fr] = readMulticolumn('Dir', Dir, 'File', File, ...
'column', 1:max(cell2mat(SET)),...
'skip', 2, 'output', 'normalized+1');
end
% select only the aminoacids you want to see
cAA = char(AA);
AA=cellstr(cAA(:,SHOWaminoACIDS));
SQUARE=zeros(size(Fr,1),1);
IX=zeros(size(Fr,1),numel(CONTROL_SETS));
i=0;
for j=CONTROL_SETS
i=i+1;
ratio(:,i) = mean(Fr(:,SET{TEST_SET}), 2) ./ mean(Fr(:,SET{j}), 2);
[~,confi(:,i)] = ttest2(Fr(:,SET{TEST_SET})',Fr(:,SET{j})',....
p_cutoff,'both','unequal');
IX(:,i) = (confi(:,i) <= p_cutoff) & (ratio(:,i) >= R_cutoff);
SQUARE = SQUARE + ratio(:,i).^2;
if PLOT_VOLCANO
subplot(1,numel(CONTROL_SETS),i);
plot(log2 (ratio(:,i)),...
-log10(confi(:,i)),'d',...
'MarkerSize',4,...
'MarkerFaceColor',0.5*[1 1 1],...
'MarkerEdgeColor',0.5*[1 1 1]); hold on;
plot( log2 (ratio(find(IX(:,i)),i)),...
-log10(confi(find(IX(:,i)),i)),'d',...
'MarkerSize',4,...
'MarkerFaceColor','r',...
'MarkerEdgeColor','r'); hold on;
line([log2(R_cutoff) MaxX],[-log10(p_cutoff) -log10(p_cutoff)]);
line([ log2(R_cutoff) log2(R_cutoff)],...
[-log10(p_cutoff) -log10(vert_cutoff)]);
xlim([-MaxX MaxX]);
end
end
R2 = sqrt(SQUARE);
IXall = find((sum(IX,2)==size(IX,2))); %hits that satisfy all criteria
%IXall = find( (sum(IX,2)>=5) ); %hits that satisfy 5 criteria
hits = char(AA(IXall,:));
Rhits = ratio(IXall,:);
R2hits = R2(IXall);
%%%%%%%%%%%% this is part where hits are clustered by H-dist %%%%%%%%%%%%%%
% figure(2)
% if size(hits,1)>3
% Y = pdist(hits,'hamming');
% Z = linkage(Y,'complete');
% [H,T,perm] = dendrogram(Z,0,'colorthreshold',20);
% set(H,'LineWidth',2)
% for i =1:size(hits,1)
% label{i} = i;
% end
% set(gca,'XTick', 1:1:size(hits,1), 'XTickLabel',label);
% hits = hits(perm,:);
% Rhits = Rhits(perm,:);
% R2hits = R2(perm);
% IXall = IXall(perm);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display all results as heat map
figure(3)
if size(IXall,1)>=HITS2DISPLAY
N=HITS2DISPLAY; % display only the first or defined number of hits
else
N=size(IXall,1); %display all
end
FrPPM = round(10^6*Fr); % convert normalized fraction frequency to PPM
imagesc( log10([FrPPM(IXall(1:N),:) ratio(IXall(1:N),:) ]+1) );
set(gca,'YTick', 1:1:N, 'YTickLabel',cellstr(hits(1:N,:)),'TickDir','out',...
'FontName','Courier New','FontSize',14);
set(gca,'XTick', 1:1:size(Fr,2)+4, 'TickDir','out');
jet1=jet;
jet1(1,:)=[0.4 0.4 0.4];
colormap(jet1);
colorbar;
% generate a plain text table for saving or copy from command line
S = char(32*ones(size(hits,1),2));
L = [ S(:,1) char(124*ones(size(hits,1),1)) S(:,1)];
F = FrPPM(IXall,:); % display frequency in ppm
toSave = [hits S ];
for i=1:numel(SET)
for j=1:numel(SET{i})
toSave = [toSave num2str(F(:,SET{i}(j))) S];
end
toSave = [toSave L];
end
toSave = [toSave S num2str(round(Rhits)) L];
disp(toSave);
if ~isempty(toSave)
fs = fopen(fullfile(Dir,'DEanalysis'),'w');
RET = char(10*ones(size(toSave,1),1));
fprintf( fs,'%s\r\n', [toSave RET]');
fclose all;
end
% or you can just copy paste the results from the command line