function FDwave_model_derived_3Delastic_2N(varargin)

% This function calculates elastic (lame's) parameters
% Complete Syntax:
%       model_derived_elastic_g1('VP',value, 'VS',value, 'RHO',value)
% Description of parameters:
%       VP           : Velocity of P-wave
%       VS           : Velocity of S-wave
%       RHO          : Density
%
%       Note: All parameters are mandatory.
%       In case any parameter is not provided the program will try to use the values stored in input folder.
%       If it cannot find any stored value then it shows error.
% Example:
%       model_derived_elastic_g1('vp',2300,'vs',2100,'rho',1600)
%       model_derived_elastic_g1();
% The assumed grid arrangement stress and velocity is :
%                   txx,tzz----------vx----------txx,tzz------>
%                  lbd,mu  |         bh             |
%                          |          .             |
%                          |          .             |
%                          |          .             |
%                      vz  |..........txz           |
%                      bv  |          muvh          |
%                          |                        |
%                          |                        |
%                          |                        |
%                   txx,tzz|------------------------|
%                          |
%                         \|/
%                          |


for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';        wfp=varargin{i+1};
        case 'vp';                  vpm=varargin{i+1};
        case 'vs';                  vsm=varargin{i+1};
        case 'rho';                 rhom=varargin{i+1};
        case 'verbose';             verbose=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';
str1=str0;          str2=str0;      str3=str0;        str4=str0;      str5=str0;    str6=str0;
if ~exist('wfp','var');        wfp=pwd;      end
if ~exist('vpm','var');        load([wfp,[filesep,'model']],'vpm');    str1=[str0,'(Default, stored)'];      end;
if ~exist('vsm','var');        load([wfp,[filesep,'model']],'vsm');    str2=[str0,'(Default, stored)'];      end;
if ~exist('rhom','var');       load([wfp,[filesep,'model']],'rhom');   str3=[str0,'(Default, stored)'];       end;
if ~exist('C44','var');        load([wfp,[filesep,'model']],'C44');    str4=[str0,'(Default, stored)'];      end;
if ~exist('C55','var');        load([wfp,[filesep,'model']],'C55');    str5=[str0,'(Default, stored)'];      end;
if ~exist('C66','var');        load([wfp,[filesep,'model']],'C66');   str6=[str0,'(Default, stored)'];       end;
if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
    disp('    FUNC: Deriving Elastic Parameter for grid type g1');
    disp([str0,' Vp model supplied',str1] )
    disp([str0,' Vs model supplied',str2] )
    disp([str0,' Density model supplied',str3] )
    disp([str0,' C44 supplied',str4] )
    disp([str0,' C55 supplied',str5] )
    disp([str0,' C66 supplied',str6] )
end

mu = vsm.^2.*rhom;
lamu = vpm.^2.*rhom;
lam = lamu-2*mu;

b=1./rhom;
bx=.5*(b(:,1:end-1,:)+b(:,2:end,:));
by=.5*(b(:,:,1:end-1)+b(:,:,2:end));
bz= .5*(b(1:end-1,:,:)+b(2:end,:,:));


muxyz = 0.125*( mu(2:end-1,2:end-1,2:end-1) + mu(3:end,2:end-1,2:end-1) + mu(2:end-1,3:end,2:end-1) +...
                mu(2:end-1,2:end-1,3:end) +mu(3:end,3:end,2:end-1) + mu(3:end,2:end-1,3:end) + mu(2:end-1,3:end,3:end) +...
                mu(3:end,3:end,3:end) );
C44 = 0.125*( C44(2:end-1,2:end-1,2:end-1) + C44(3:end,2:end-1,2:end-1) + C44(2:end-1,3:end,2:end-1) +...
                C44(2:end-1,2:end-1,3:end) +C44(3:end,3:end,2:end-1) + C44(3:end,2:end-1,3:end) + C44(2:end-1,3:end,3:end) +...
                C44(3:end,3:end,3:end) );
C55 = 0.125*( C55(2:end-1,2:end-1,2:end-1) + C55(3:end,2:end-1,2:end-1) + C55(2:end-1,3:end,2:end-1) +...
                C55(2:end-1,2:end-1,3:end) +C55(3:end,3:end,2:end-1) + C55(3:end,2:end-1,3:end) + C55(2:end-1,3:end,3:end) +...
                C55(3:end,3:end,3:end) );
C66 = 0.125*( C66(2:end-1,2:end-1,2:end-1) + C66(3:end,2:end-1,2:end-1) + C66(2:end-1,3:end,2:end-1) +...
                C66(2:end-1,2:end-1,3:end) +C66(3:end,3:end,2:end-1) + C66(3:end,2:end-1,3:end) + C66(2:end-1,3:end,3:end) +...
                C66(3:end,3:end,3:end) );   
				
str=strcat(wfp,[filesep,'derived_param']);
save(str,'lam','lamu','mu','b','bx','by','bz','muxyz','C44','C55','C66')

if strcmp(verbose,'y')
    disp([str0,' Derived parameters generated successfully'])
    disp([str0,' Derived parameters are saved in "',str])
end
end

