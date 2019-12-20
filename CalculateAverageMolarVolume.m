load('CustomDataForLatticeModelBovine.mat');
whichIs=24;
sequence=protein(whichIs).tail.Sequence;
load('AminoAcidVolumes.mat')
AA='E';
AA_vol=ones(length(sequence),1);
for j=1:length(sequence)
%    strfind
%    sequence(j) 
    AA=upper(sequence(j));
    ind=cellfun(@(x)strcmp(AA,x),Letter);
    AA_vol(j)=Partial_volume(ind);
end
Avg=sum(AA_vol)/length(AA_vol)
bar(AA_vol)