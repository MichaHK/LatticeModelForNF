function []=main_bulk(handles)

1;
temp=get(handles.text5,'String');
FileList=get(handles.listbox1,'String');
filename=FileList{get(handles.listbox1,'Value')};
axes(handles.axes1);
1;


Exp_Pressure=importdata(filename);
DistanceData=Exp_Pressure(:,1)/10;
PressureData=Exp_Pressure(:,2);
1;
scatter(DistanceData,PressureData);
set(gca,'yscale','log');

xlim([10,90])
ylim([10,10^7])


xx1 = smooth(PressureData,DistanceData,0.01,'rloess');
% pp=spline(PressureData,DistanceData);
pp=csaps(DistanceData,PressureData,1);
pp2=csaps(PressureData,DistanceData,0.2);
pp3=csaps(log10(PressureData),DistanceData,0.99);

PEG_Data=abs((log10((PressureData-130).*10^-0.57)/2.75).^(1/0.21));
pp_Dist_PEG=csaps(PEG_Data,DistanceData,0.8); % @D(wt)
p_der_PEG=fnder(pp_Dist_PEG,1); % @D'(wt)
NewPegData=linspace(PEG_Data(1),PEG_Data(end),20); % wt: new wt values for interpolation of distance D
NewX=linspace(DistanceData(1),DistanceData(end-1));
Pressure_PEG=10.^(0.57+2.75*NewPegData.^0.21)+130; %P(wt)
Deriv_Pressure_PEG=log(10).*10.^(0.57+2.75*NewPegData.^0.21).*(2.75*0.21*NewPegData.^-1.21);
Bulk_wt=-ppval(pp_Dist_PEG,NewPegData)/2.*...
        ppval(p_der_PEG,NewPegData).*...
        Deriv_Pressure_PEG;



% use this with new peg to plot. 
    
    %     a-\frac{b}{D}+\left[\frac{1}{cD^{2}+d}\right]^{2}

p_der=fnder(pp,1);
% NewX=DistanceData;
NewX=linspace(DistanceData(1),DistanceData(end-1));
NewY2=linspace(PressureData(1),PressureData(end));
1;
y_prime=ppval(p_der,NewX);
y_prime_PEG=ppval(p_der_PEG,NewY2);
PressureData_PEG=y_prime_PEG;
% Bulk=-NewX/2.*y_prime;
% Bulk_PEG=-NewX/2.*

hold on;
plot(NewX,ppval(pp,NewX))
plot(ppval(pp_Dist_PEG,NewPegData),Pressure_PEG,'r') % D(wt). For interpolated D on a new set of wt. 

% plot(ppval(pp3,log10(NewY2)),NewY2,'r')
% plot(xx1,PressureData,'r')
% plot(ppval(pp2,NewY2),NewY2,'r')
hold off;


axes(handles.axes2);
% scatter(NewX,Bulk);
% set(gca,'yscale','log');

1;

% ExecuteLine=get(handles.Command_edit,'String');
% eval(ExecuteLine);




    function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function out=SetUnits(Pressure) %from nm^-3 to Pa
        Pressure=Pressure*10^27;
        Kb= 1.38065*10^-23; % Boltzmann constant's: J/K  (joules per kelvin) 
        T= 298.15; % Temperature: Kelvins. 
        out=Pressure*Kb*T; 
    end


end



