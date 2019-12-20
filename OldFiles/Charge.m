function [AACharges,TotalCharge,PartialCharge,PartialNegative, PartialPositive, ListNegSites, ListPosSites]=...
    Charge(Sequence,PhosIndices,pH)
% phos indices includes all indices (phospho*...) 
% PhosIndices are indexed according to tail (region) indices, not the
% complete protein indices. 


pKset = {'N_term',8.6;'K',10.8;'R',12.5;'H',6.5;'D',3.9;'E',4.1;...
     'C',8.5;'Y',10.1;'C_term',3.6};
Sequence=upper(Sequence);

PhosCharge=0; 
TermChargeCorr=0;
IncludeTerminals=0;   

    if nargin<2;
        PhosIndices=[];
        pH=7;
    end
 
    AACharges=zeros(1,length(Sequence));
    for i=1:numel(AACharges)
        AACharges(i)=AACharge(Sequence(i),pH,pKset);
    end
    
    1;
    
    if ~isempty(PhosIndices)
        AACharges(PhosIndices)=AACharges(PhosIndices)-2; 
    end

TotalCharge=sum(AACharges);
PartialCharge=TotalCharge/numel(Sequence);

[PartialNegative, PartialPositive, ListNegSites, ListPosSites]= CalcPartNegPos(AACharges,Sequence);

% plot(AACharge,'b');
% hold on;
% plot(smooth(AACharge,5),'r');

end


function Charge=AACharge(L,pH,pKset)
%     pKset = {'N_term',8.6;'K',10.8;'R',12.5;'H',6.5;'D',3.9;'E',4.1;...
%             'C',8.5;'Y',10.1;'C_term',3.6};
    Charge=0;
    acids='CDEY';  %C_term also. 
    bases='KRH';
    if ~isempty(strfind(acids,L))
        LetterIndInpKset=find(cellfun(@(c) ~isempty(strfind(c,L)),pKset(:,1)));
        pKa=pKset{LetterIndInpKset,2};
        FracCharge=1/(1+10^(pKa-pH));
        Charge=-FracCharge;
    end
    if ~isempty(strfind(bases,L))
        LetterIndInpKset=find(cellfun(@(c) ~isempty(strfind(c,L)),pKset(:,1)));
        pKa=pKset{LetterIndInpKset,2};
        FracCharge=1/(1+10^(-pKa+pH));
        Charge=FracCharge;
    end
end

function [PartialNeg, PartialPos, ListNegSites, ListPosSites]= CalcPartNegPos(AACharges,Sequence)
    ListNegSites=find(AACharges<0);
    ListPosSites=find(AACharges>0);
    PartialNeg=sum(AACharges(ListNegSites))/length(Sequence);
    PartialPos=sum(AACharges(ListPosSites))/length(Sequence);
end
