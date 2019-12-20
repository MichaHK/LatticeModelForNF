function ChargeDistNFLtip
if (~exist('NFLtipComparison.mat','file'))
    ProteinsChart=importdata('NFLtipComparison.xlsx');
    save('NFLtipComparison.mat','ProteinsChart');
else
    load('NFLtipComparison.mat')
end
Chart=ProteinsChart.textdata.Sheet1;
[Length,~]=size(Chart);
PhosIndices=[];
AverageWindow=6;
FigH=figure;
BovineSeq='MSSFSYEPYYSTSYKRRYVETPRVHISSVRSGYSTARSAYSSYSAPVSSSLSVRRSYSSSSGSLMPSLESLDLSQVAAISNDLKSIRTQEKAQLQDLNDRFASFIERVHELEQQNKVLEAELLVLRQKHSEPSRFRALYEQEIRDLRLAAEDATNEKQALQGEREGLEETLRNLQARYEEEVLSREDAEGRLMEARKGADEAALARAELEKRIDSLMDEIAFLKKVHEEEIAELQAQIQYAQISVEMDVSSKPDLSAALKDIRAQYEKLAAKNMQNAEEWFKSRFTVLTESAAKNTDAVRAAKDEVSESRRLLKAKTLEIEACRGMNEALEKQLQELEDKQNADISAMQDTINKLENELRTTKSEMARYLKEYQDLLNVKMALDIEIAAYRKLLEGEETRLSFTSVGSLTTGYTQSSQVFGRSAYGGLQTSSYLMSARSFPSYYTSHVQEEQIEVEETIEAAKAEEAKDEPPSEGEAEEEEKEKEEAEAEAEAEAEAEAEEEEGAQEEEAAKEDAEEAKEEEGGEGEEAEETKEAEEEEKKDEGAGEEQATKKKD';
BocIndices=(numel(BovineSeq)-110):numel(BovineSeq);
for i=2:Length
    Sequence=Chart{i,5};
    Name=Chart{i,3};
    Class=Chart{i,6};
%     Identities=Chart{i,7};
%     Positives=Chart{i,8};

    SeqLength=numel(Sequence);
    [AACharges,TotalCharge,PartialCharge,PartialNegative, PartialPositive, ListNegSites, ListPosSites]=ChargeCalculate(Sequence,PhosIndices);
    indices=(numel(Sequence)-110):numel(Sequence);

    AACharges=smooth(AACharges,AverageWindow);
    [ListNegSites ListPosSites]=ListNegPos(AACharges);
    subplot(2,5,i-1)
    hold on;
    set(gca,'fontsize',14)
    % ylabel ({'Charge' 'density'});
    bar(ListNegSites,AACharges(ListNegSites),'FaceColor','r','EdgeColor','r')
    bar(ListPosSites,AACharges(ListPosSites),'FaceColor','b','EdgeColor','b')
    % set(gca,'Linewidth',1.5);
    AAindex=hydrophobicity(Sequence,0,1);
    scatter(find(AAindex>1),zeros(numel(find(AAindex>1)),1),20,'*','k'); % FIX TO INCLUDE STARS FOR HYDROPHOBICITY
    [GlobalScore, GlobalAlignment]=nwalign(BovineSeq(BocIndices),Sequence(indices));
%     title({Name,[Class ' , ' num2str(SeqLength) 'aa, (' Identities ',' Positives ')%']},'FontSize',12)
    title({Name,[Class ' , ' num2str(SeqLength) 'aa, ' num2str(GlobalScore,3)]},'FontSize',12)
%     title([ Name ', ' num2str(SeqLength) 'aa'],'FontSize',12)
%     showalignment(GlobalAlignment)
    xlim([numel(Sequence)-110 numel(Sequence)])
    set(gca,'FontSize',12)
    ylim([-1 0.5])
    
end

% AxesH    = findobj(FigH, 'Type', 'Axes');
% YLabelHC = get(AxesH, 'YLabel');
% YLabelH  = [YLabelHC{:}];
% set(YLabelH, 'String', 'Y-label')
% TitleHC  = get(AxesH, 'Title');
% TitleH   = [TitleHC{:}];
% set(TitleH, 'String', 'The title');

function [ListNegSites ListPosSites]= ListNegPos(AACharge)
ListNegSites=find(AACharge<0);
ListPosSites=find(AACharge>0);
end
end