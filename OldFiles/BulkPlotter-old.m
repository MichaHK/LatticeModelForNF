% Fit, Bulk
ALLFiles={'L-deH_150mM.csv',  'LH_150mM_MAX2015.csv',  'L-deM_150mM.csv', 'LM_150mM_MAX2015.csv','de-LMH_150mM.csv','LMH_150mM_MAX2015.csv'}; % 'LMH_150mM_Roy.csv','LH_150mM.csv',
LMH_Files={'de_LMH_150mM.csv','LMH_150mM_MAX2015.csv'}; % 'LMH_150mM_Roy.csv'
LH_Files= {'L-deH_150mM.csv',  'LH_150mM_MAX2015.csv'}; % 'LH_150mM.csv'
LM_Files= {'L-deM_150mM.csv', 'LM_150mM_MAX2015.csv'}; % 'LM_150mM.csv',
Phos_files={'LMH_150mM_MAX2015.csv','LH_150mM_MAX2015.csv','LM_150mM_MAX2015.csv','L_150mM.csv'};% 'LMH_150mM_Roy.csv'
DePhos_files={'L-deH_150mM.csv','L_deM_150mM.csv','de_LMH_150mM.csv'};
% files=LMH_Files;
% SmoothingParams=[0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95];
 
% files=LH_Files;
% SmoothingParams=[0.95,0.95,0.005,0.95,0.95,0.95,0.95,0.95];
 
% files=LM_Files;
% SmoothingParams=[0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95];

% files=Phos_files;
% SmoothingParams=[0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95];

files=DePhos_files;
SmoothingParams=[0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95];

% files=ALLFiles;
% SmoothingParams=[0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95];



BulkArray=[];
NewDistancesArray=[];
figure( 'Name', 'Fits' );
for i=1:numel(files)
%     for i=1:2
        Data=importdata(files{i});
        Pressure=Data(:,2);
        Distance=Data(:,1);
        h=subplot(numel(files),2,2*i-1);
        [NewDistances,NewPressures,bulk] = createFitDefault(Pressure,Distance,SmoothingParams(i));  
        BulkArray=[BulkArray bulk];
        NewDistancesArray=[NewDistancesArray NewDistances];
        s=semilogy(NewDistances/10,NewPressures);
        hold on; 
        scatter(Distance/10,Pressure,'r','o');
        xlim([15 100])
        ylim([60 10^7])
        % Label axes
         if (i==numel(files)) xlabel( 'Interfilament distance [nm]' ); end
        ylabel( 'Pressure [Pa]' );
        set(gca,'Ytick',[10.^(2:2:7)])
        set(gca,'Xtick',[20:20:100])

        grid on;       
        DisplayName=strrep(files{i},'_MAX2015','');
        DisplayName=strrep(DisplayName,'_150mM','');
        DisplayName=strrep(DisplayName,'.csv','');
        DisplayName=strrep(DisplayName,'_',' ');
        legend( h, ['Fit for ' DisplayName] , 'Experiment', 'Location', 'NorthEast' );
        g=subplot(numel(files),2,2*i);
        semilogy(NewDistances(2:end)/10,bulk)
        legend( g, 'Bulk modulus' , 'Location', 'NorthEast' );
        xlim([15 100])
        ylim([130 10^8])
        set(gca,'Ytick',[10.^(2:2:7)])
        set(gca,'Xtick',[20:20:100])
        % Label axes
        if (i==numel(files)) xlabel( 'Interfilament distance [nm]' ); end
        ylabel( 'Osmotic bulk pressure [Pa]' );
        grid on;

end

figure ('Name', 'Bulks');
colors={'b','r','k','g','m'};
for i=1:numel(files)
    NewDistances=NewDistancesArray(:,i);
    semilogy(NewDistances(1:end-1)/10,BulkArray(:,i),colors{i})
    hold on;
end
1;
% Label axes
        xlabel( 'Distance [nm]' );
        ylabel( 'Osmotic Bulk Mod. [Pa]' );
        xlim([15 100])
        ylim([130 10^9])
        grid on;
legend(gca, files , 'Location', 'NorthEast' );