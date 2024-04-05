function temp2 = cell2numaa(cNuc4, Fr4, AA)
% this functiona converts a string of peptides and their frequencies to
% a 400x400 array of numbers. Peptides are in "cNuc4"-variable.

    % convert amino acids to numbers that reflect their order in AA string
    nNuc4=zeros(size(cNuc4));
    NN=numel(AA);

    for i=1:NN
        IX=find(cNuc4==AA(i));
        nNuc4(IX)=i;
    end

    % eliminate peptides that have unnatural characters
    [row,~]=find(nNuc4==0);
    IX2 = setdiff(1:size(nNuc4,1),row);
    nNuc4 = nNuc4 ( IX2, :);
    Fr4 = Fr4 ( IX2, :);
    
    %X position contains 1st and 3rd letter, Y position - 2nd and 4th
    X = NN*(nNuc4(:,1)-1) + nNuc4(:,3);
    Y = NN*(nNuc4(:,2)-1) + nNuc4(:,4);

    temp2=zeros(NN^2,NN^2);


    for i=1:numel(Fr4);
        temp2( X(i), Y(i) ) = Fr4(i);
    end
    
end