function wiggle (varargin)
%WIGGLE plots seismic data using wiggles, and the willges and positive lobes could be
% one of the three different colors, i.e., black ('k'), blue('b'), and red ('r').

% WIGGLE(X) plots each column of X as curves vertically. 
% WIGGLE(X,WC) plots each column of X as curves vertically, the wiggle color is defined by WC. 
% WIGGLE(X,WC,LC) plots each column of X as curves vertically, 
%                 the wiggle color is defined by WC and the lobe color is defined by LC. 


%   Wiggles Color: 'k','b','r'
%     the default is black ('k'). 
%
%   Lobes color: 'k','b','r'
%     the default is black ('k'). 

switch (nargin)
    
    case 0
        error('No input arguments!');
                
        % Data and properties
    case 1
        data = check_data(varargin{1});
        [nz,nx] = size(data);
        x = 1:nx;
        z = 1:nz;
        wc=str2rgb('k');
        lc=str2rgb('k');
    case 2
        data = check_data(varargin{1});
        [nz,nx] = size(data);
        x = 1:nx;
        z = 1:nz;
        if ischar(varargin{2})
            wc=str2rgb(varargin{2});
        else
            error('Wiggle: properties must be a string.');
        end
        lc=str2rgb('k');
    case 3
        data = check_data(varargin{1});
        [nz,nx] = size(data);
        x = 1:nx;
        z = 1:nz;
        if ischar(varargin{2}) && ischar(varargin{3})
            wc=str2rgb(varargin{2});
            lc=str2rgb(varargin{3});
        else
            error('Wiggle: color property must be a string.');
        end
        
    otherwise
        error('Too many input arguments!');
        
end
   
trmx= max(abs(data));
amx=mean(trmx);

% take the average as dx
dx1 = abs(x(2:nx)-x(1:nx-1));
dx = median(dx1);

dz=z(2)-z(1);
data= data * dx /amx; 

% set display range 
x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
z1=min(z)-dz; z2=max(z)+dz;
 
set(gca,'NextPlot','add','Box','on','XLim', [x1 x2],...
    'YDir','reverse','YLim',[z1 z2],'xAxisLocation','top');
 
linewidth = 0.1;
z=z'; 	% input as row vector
zstart=z(1);
zend  =z(nz);

for i=nx:-1:1
   
	tr=data(:,i); 	% one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% zero crossing points
	%npos = length(i1);


	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); % locations with zero amplitudes
	aadd = zeros(size(zadd));

	[zpos,~] = find(tr >0);
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
	if tr(1)>0 	
        a0=0; 
        z0=1.00;
    else
        a0=0; 
        z0=zadd(1);
	end
    
    if tr(nz)>0	
        a1=0; 
        z1=nz;
    else
        a1=0; 
        z1=max(zadd);
    end
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		

	zzz = zstart + zz*dz -dz;
    patch( aa+x(i) , zzz,  lc);
    line( 'Color',[1 1 1],'Xdata', x(i)+[0 0], 'Ydata',[zstart zend]); % remove zero line
    for i=1:nx
      if trmx(i) ~= 0    % skip the zero traces
	     tr=data(:,i); 	
	     line( 'Color',wc,'LineWidth',linewidth, ...
	     'Xdata', tr+x(i), 'Ydata',z);	% negatives line
      end
    end
end
 
 
 %%
function data = check_data(d)

if isnumeric(d)
    data = d;
else
    error('Wiggle: data must be numeric');
end

if size(d,1)<2||size(d,2)<2
    error('Wiggle: data size must be larger than 1');
end

end

%%
function rgb = str2rgb(s)
switch s    
    case 'k'
        rgb = [0, 0, 0];
    case 'b'
        rgb = [0, 0, 1];
    case 'r'
        rgb = [1, 0, 0];
    otherwise
        rgb = NaN;
end
end

end

