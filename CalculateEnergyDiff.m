function CalculateEnergyDiff
% Do energy per monomer or phos site!!!!!
LMH_Files={'de-LMH_150mM.csv','LMH_150mM_MAX2015.csv'}; % 'LMH_150mM_Roy.csv'
LH_Files= {'L-deH_150mM.csv',  'LH_150mM_MAX2015.csv'}; % 'LH_150mM.csv'
LM_Files= {'L-deM_150mM.csv', 'LM_150mM_MAX2015.csv'}; % 'LM_150mM.csv',
type='LM'; 

switch type
    case 'LM'
        ExperimentalFileNames=LM_Files;
        AverageNumOfMonomersInTail=0.7*158+0.3*514; % 514 for nfm, 679 nfh
        AverageNumOfPhosSitesInTail=0.7*1+0.3*27;
    case 'LH'
        ExperimentalFileNames=LH_Files;
        AverageNumOfMonomersInTail=0.8*158+0.2*679; % 514 for nfm, 679 nfh
        AverageNumOfPhosSitesInTail=0.8*1+0.2*27;
    case 'LMH'
        ExperimentalFileNames=LMH_Files;
        AverageNumOfMonomersInTail=0.6*158+0.27*514+0.13*679; % 514 for nfm, 679 nfh
        AverageNumOfPhosSitesInTail=0.6*1+0.27*27+0.13*48;
end



% close ('all')

RodLength=450; % In Angstorms!!!!!
RodCyl=50;      % In Angstorms!!!!!
NumOfTailsInVol=16;
NumOfMonomersInVolume=AverageNumOfMonomersInTail*NumOfTailsInVol;
NumOfPhosSitesInVolume=NumOfMonomersInVolume*AverageNumOfPhosSitesInTail/AverageNumOfMonomersInTail;
FitFiles=cellfun(@(x) strcat('Fit_',x), ExperimentalFileNames, 'UniformOutput', false);
    Data1=importdata(FitFiles{1});
    Pressure1=Data1(:,2);
    Distance1=Data1(:,1);
    Data2=importdata(FitFiles{2});
    DisplayName1=FixName(FitFiles{1});
    DisplayName2=FixName(FitFiles{2});
    Pressure2=Data2(:,2);
    Distance2=Data2(:,1);
    NewDistances=linspace(max(min(Distance1),min(Distance2)),...
                  max(max(Distance1),max(Distance2)),600);
    
    NewPressure1=interp1(Distance1,Pressure1,NewDistances);
    NewPressure2=interp1(Distance2,Pressure2,NewDistances);
    NewPressure1=ReplaceNaN(NewPressure1);
    NewPressure2=ReplaceNaN(NewPressure2);
    
    DeltaPressures=NewPressure2-NewPressure1;
    %% Plot subtracted pressures
%     h=figure;
%     plot(NewDistances,DeltaPressures,'k')
% %     plot(NewDistances,DeltaPressures,'k')
%     hold on;
%     plot(NewDistances,NewPressure1,'b')
%     plot(NewDistances,NewPressure2,'r')
%     
%     
%     legend(gca,[DisplayName2 '-' DisplayName1],DisplayName1,DisplayName2)  
    %% integral Volume
    1;
%     Volume=RodLength*pi.*(NewDistances/2-RodCyl).^2;               % water(!) volume
    Volume=RodLength*(sqrt(3)*NewDistances.^2)/4-(pi*RodCyl^2)/2;
    integrad=DeltaPressures(2:end).*diff(Volume);
    TriangMat=triu(ones(numel(integrad),numel(integrad)));
    integral=TriangMat*integrad';
%     integral=cumsum(integrad(end:-1:1));
%     integral=integral(end:-1:1);
    figure
    KbT_integral=integral*(2.89719*10^-9); % http://www.wolframalpha.com/input/?i=pascal*%28angstrom%5E3%29+%2F%28Boltzmann+constant*25+Celsius%29
    if (0)
    plot(Volume(2:end)/NumOfMonomersInVolume,KbT_integral/NumOfMonomersInVolume)
    legend(gca,'Integral')
    title([DisplayName2(4:end) ': \Delta{G}(V)=\int_{V_{max}}^{V} (\Pi_{phos}-\Pi_{de})dV'])
    xlabel('Volume (A^3/monomer): space in A^3 for each monomer')
    ylabel('\Pi-V work KT/monomer')
    end
    if(1)
    V_a=0.134;  % In Angstorms^3!!!!!
    VolumeFraction=NumOfMonomersInVolume*V_a./Volume;
%     plot(VolumeFraction(2:end),KbT_integral/NumOfPhosSitesInVolume)
    plot(NewDistances(2:end)/10,KbT_integral/NumOfPhosSitesInVolume)

    legend(gca,'Integral')
    title([DisplayName2(4:end) ': \Delta{G}(V)=\int_{V_{max}}^{V} (\Pi_{phos}-\Pi_{de})dV'])
%     xlabel('\phi^{-1}')
    xlabel('d (nm)')
    ylabel('\Pi-V work KT/phosphosite')
    end
    
    %% Save Data to files
    if (0)
%     filename = [pwd '\WorkCalc\' DisplayName2 '_V_perMono_vs_E_KT_perMono' '.csv'];
%     csvwrite(filename, [Volume(2:end)/NumOfMonomersInVolume;
%     KbT_integral'/NumOfMonomersInVolume]');
    end
    if (1)
    filename = [pwd '\WorkCalc\' DisplayName2 '_Dist_vs_E_KT_perPhosSite' '.csv'];
    csvwrite(filename, [NewDistances(2:end); KbT_integral'/NumOfPhosSitesInVolume]');

    end
    %% integral distance
%     figure
%     integrad=DeltaPressures(2:end).*diff(NewDistances);
%     TriangMat=triu(ones(numel(integrad),numel(integrad)));
%     integral=TriangMat*integrad';
%     plot(NewDistances(2:end),integral)
%     legend(gca,'Integral')
%     title([DisplayName2(4:end) ': \Delta{G}(D)=\int_{D_{max}}^{D} (\Pi_{phos}-\Pi_{de})dD'])
%     xlabel('Distance[A]')
%     ylabel('\Pi-D [A\cdot{Pascal}]')
%% Aux function
    function NewPressure=ReplaceNaN(OldPressure)
        NewPressure=OldPressure;
        indices=isnan(OldPressure);
        NewPressure(indices)=0; 
    end
    
    function DisplayName=FixName(filename)
        DisplayName=strrep(filename,'_MAX2015','');
        DisplayName=strrep(DisplayName,'_150mM','');
        DisplayName=strrep(DisplayName,'.csv','');
        DisplayName=strrep(DisplayName,'_',' ');
    end
end