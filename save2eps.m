%SAVE2EPS Saves a figure as a properly cropped eps
%
%   save2eps(fileName,handle,dpi)
%
%   - fileName: Destination to write the eps to.
%   - handle: (optional) Handle of the figure to write to a eps.  If
%             omitted, the current figure is used.  Note that handles
%             are typically the figure number.
%   - dpi: (optional) Integer value of dots per inch (DPI).  Sets
%          resolution of output eps.
%
%   Saves figure as a eps with margins cropped to match the figure size.
function save2eps(fileName,handle,dpi)

% Verify correct number of arguments
error(nargchk(0,3,nargin));

% Set omitted arguments
if nargin < 3
    dpi = 300;
    if nargin < 2
        handle = gcf;
        if nargin < 1
            [fileName,pathName] = uiputfile('*.eps','Save to EPS file:');
            if fileName == 0; return; end
            fileName = [pathName,fileName];
        end
    end            
end

% Backup previous settings
prePaperType = get(handle,'papertype');
prePaperUnits = get(handle,'paperunits');
preUnits = get(handle,'units');
prePaperSize = get(handle,'papersize');

% Make changing paper type possible
set(handle,'papertype','<custom>');

% Set units to all be the same
set(handle,'paperunits','inches');
set(handle,'units','inches');

% Set the page size and position to match the figure's dimensions
position = get(handle,'position');
set(handle,'paperposition',[0,0,position(3:4)]);
set(handle,'papersize',position(3:4));

% Save the eps (this is the same method used by "saveas")
print(handle,'-depsc2',fileName,sprintf('-r%d',dpi))

% Restore the previous settings
set(handle,'papertype',prePaperType);
set(handle,'paperunits',prePaperUnits);
set(handle,'units',preUnits);
set(handle,'papersize',prePaperSize);