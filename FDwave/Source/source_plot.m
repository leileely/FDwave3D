function source_plot(varargin)%figno,spind1,spind2,spind3)

% This function generates the Sine wavelet signature.
% Complete Syntax:
%       source_plot('WFP',path,'Wave_Type',option,'FigNO',value,'I1',value,'I2',value,'I3',value)
% Description of parameters:
%       WFP         : Path to working directory
%       Wave_Type   : Type of wave  e.g. Elastic 
%       FigNO       : Figure number
%       I1,I2,I3    : subfigure index
% Example:
%       source_plot('WFP',pwd,'FigNO',1,'I1',1,'I2',1,'I3',1)

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';             wfp=varargin{i+1};
        case 'wave_type';       wave_type = varargin{i+1};
        case 'figno';           figno = varargin{i+1};
        case 'i1';              spind1 = varargin{i+1};
        case 'i2';              spind2 = varargin{i+1};
        case 'i3';              spind3 = varargin{i+1};
        otherwise;              error('%s is not a valid argument name',varargin{i});
    end
end

if ~exist('wfp','var');             wfp=pwd;                end;
if ~exist('wave_type','var');       wave_type='Elastic';    end
if ~exist('figno','var');           figno=figure();         end

if ~exist('spind1','var')&&~exist('spind2','var')&&~exist('spind3','var');
     spind1=1;  spind2=1;   spind3=1;
end

if ~exist('spind1','var')||~exist('spind2','var')||~exist('spind3','var');
    if strcmp(wave_type,'Acoustic1');            spind1=1;  spind2=2;   spind3=2;
    elseif strcmp(wave_type,'Acoustic2');        spind1=2;  spind2=2;   spind3=3;
    elseif strcmp(wave_type,'Elastic');          spind1=2;  spind2=2;   spind3=4;
    elseif strcmp(wave_type,'Viscoelastic');     spind1=2;  spind2=3;   spind3=6;
    end
end

str_s=[wfp,[filesep,'source']];        load(str_s);
str_m=[wfp,[filesep,'model']];         load(str_m,'wave_type')

figure(figno);
h=subplot(spind1,spind2,spind3);


plot(src(1:tdn));

if verLessThan('matlab','8.4')
    xtk = dt*get(h,'XTick');        set(h,'XTickLabel',xtk)
    ytk = 1*get(h,'YTick');        set(h,'YTickLabel',ytk)
    ax = gca;                       set(ax,'UserData',[dt,1]);
    h=gcf;                          set(h,'ResizeFcn','axTickUpdate');
else 
    %%% works in MATLAB 2014 above
    ax = gcf;                       
    ax.UserData.dh = 1;         ax.UserData.dv = 1;
    ax.SizeChangedFcn = 'axTickUpdate';
end

title(['Source (',src_name,')']);
xlabel('Time (seconds)');   ylabel('Amplitude');



