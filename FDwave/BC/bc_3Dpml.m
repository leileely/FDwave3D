function BC=bc_3Dpml( cpath,BCtype,nAB )
% BC values of PML condition.
%   nAB : no of nodes used for damping
%   BCtype : type of BC viz. 'topFS' means top edge as free surface
%                            'allABC' means all sides as absorbing
%   default : nAB=25 & BCtype='allABC'
% switch nargin
%     case 0
%         nAB = 50;
%         BCtype='topFS';
%        
%     case 2
% %         disp(['    Input ABL layer nodes are       :  ', num2str(nAB)]);
% %         disp(['    Input type of top condition is  :  ', BCtype]);
%     otherwise
%         disp('Check input arguments')
% end


% disp(['    Default ABL layer nodes are       :  ', num2str(nAB)]);
% disp(['    Default type of top condition is  :  ', BCtype]);

%% see eq. (12) of the related manuscript
BCname='PML';

load([cpath,[filesep,'model']],'nx','ny','nz','dx','dy','dz','vpm');

i=1:nAB;
wt = exp(-(0.0053*(nAB-i)).^2);

BCx = zeros(nz,nx,ny);
BCy = zeros(nz,nx,ny);
BCz = zeros(nz,nx,ny);
R=1e-6;
vm=max(vpm(:));
ppml = -log(R)*3*vm/(2*nAB^3);

if strcmp(BCtype,'topFS')
    for k=1:length(wt)
        % no top face & top corners no  BCz & top edges no BCz
        %% corners
        BCx(1:nAB,k,1:nAB)= (nAB-k+1)^2*ppml/dx; %Top 
        BCy(1:nAB,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,1:nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,nx-nAB+1:nx,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,k,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,1:nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,nx-nAB+1:nx,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,1:nAB)= (nAB-k+1)^2*ppml/dx; %Bottom 
        BCy(nz-nAB+1:nz,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,1:nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nx-nAB+1:nx,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,1:nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nx-nAB+1:nx,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
                
        %% edges
        
        BCy(1:nAB,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,nAB+1:nx-nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCy(nz-nAB+1:nz,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nAB+1:nx-nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCy(1:nAB,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        %BCz(k,nAB+1:nx-nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCy(nz-nAB+1:nz,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nAB+1:nx-nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        
        BCx(1:nAB,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        %BCz(k,1:nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(nz-k+1,1:nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        %BCz(k,nx-nAB+1:nx,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(nz-k+1,nx-nAB+1:nx,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        
        BCx(nAB+1:nz-nAB,k,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        %% faces
        BCx(nAB+1:nz-nAB,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCx(nAB+1:nz-nAB,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        
        BCy(nAB+1:nz-nAB,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy; 
        BCy(nAB+1:nz-nAB,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        %BCz(k,nAB+1:nx-nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz; 
        BCz(nz-k+1,nAB+1:nx-nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        
    end
    
elseif strcmp(BCtype,'topABC')
    for k=1:nAB
        
        %% corners
        BCx(1:nAB,k,1:nAB)= (nAB-k+1)^2*ppml/dx; %Top 
        BCy(1:nAB,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(k,1:nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        BCz(k,nx-nAB+1:nx,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,k,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(k,1:nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(1:nAB,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(k,nx-nAB+1:nx,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,1:nAB)= (nAB-k+1)^2*ppml/dx; %Bottom 
        BCy(nz-nAB+1:nz,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,1:nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nx-nAB+1:nx,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,1:nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nz-nAB+1:nz,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nx-nAB+1:nx,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
                
        %% edges
        
        BCy(1:nAB,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(k,nAB+1:nx-nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCy(nz-nAB+1:nz,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nAB+1:nx-nAB,1:nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCy(1:nAB,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(k,nAB+1:nx-nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        BCy(nz-nAB+1:nz,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        BCz(nz-k+1,nAB+1:nx-nAB,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dz;
        
        
        BCx(1:nAB,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(k,1:nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(nz-k+1,1:nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(1:nAB,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(k,nx-nAB+1:nx,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        BCx(nz-nAB+1:nz,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCz(nz-k+1,nx-nAB+1:nx,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
        
        BCx(nAB+1:nz-nAB,k,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,1:nAB,k)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,k,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,1:nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,nx-k+1,1:nAB)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,nx-nAB+1:nx,k)= (nAB-k+1)^2*ppml/dy;
        
        BCx(nAB+1:nz-nAB,nx-k+1,ny-nAB+1:ny)= (nAB-k+1)^2*ppml/dx; 
        BCy(nAB+1:nz-nAB,nx-nAB+1:nx,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        %% faces
        BCx(nAB+1:nz-nAB,k,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        BCx(nAB+1:nz-nAB,nx-k+1,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dx; 
        
        BCy(nAB+1:nz-nAB,nAB+1:nx-nAB,k)= (nAB-k+1)^2*ppml/dy; 
        BCy(nAB+1:nz-nAB,nAB+1:nx-nAB,ny-k+1)= (nAB-k+1)^2*ppml/dy;
        
        BCz(k,nAB+1:nx-nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz; 
        BCz(nz-k+1,nAB+1:nx-nAB,nAB+1:ny-nAB)= (nAB-k+1)^2*ppml/dz;
        
    end
    
end


BC=struct('pmlx',BCx,'pmly',BCy,'pmlz',BCz);
end

