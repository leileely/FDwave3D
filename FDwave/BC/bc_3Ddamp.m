function BC=bc_3Ddamp( cpath,BCtype,nAB )
% BC values of damping condition.
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

BCname='ABL';

load([cpath,[filesep,'model']],'nx','ny','nz');
BCabl = ones(nz,nx,ny);

i=1:nAB;
wt = exp(-(0.0053*(nAB-i)).^2);

if strcmp(BCtype,'topFS')
    for k=1:length(wt)
        BCabl(1:nz-k+1,k,k:ny-k+1) = wt(k); %left BC
        BCabl(nz-k+1,k:nx-k+1,k:ny-k+1)= wt(k); %Bottom BC
        BCabl(1:nz-k+1, nx-k+1,k:ny-k+1)= wt(k);    % Right BC
        BCabl(1:nz-k+1,k:nx-k+1,k)=wt(k);    %Front BC
        BCabl(1:nz-k+1,k:nx-k+1,ny-k+1)=wt(k);   %Back BC
    end
    
elseif strcmp(BCtype,'topABC')
    for k=1:length(wt)
        BCabl(k,k:nx-k+1,k:ny-k+1)= wt(k); %Top BC
        BCabl(k:nz-k+1,k,k:ny-k+1) = wt(k); %left BC
        BCabl(nz-k+1,k:nx-k+1,k:ny-k+1)= wt(k); %Bottom BC
        BCabl(k:nz-k+1, nx-k+1,k:ny-k+1)= wt(k);    % Right BC
        BCabl(k:nz-k+1,k:nx-k+1,k)=wt(k);    %Front BC
        BCabl(k:nz-k+1,k:nx-k+1,ny-k+1)=wt(k);   %Back BC
    end
    
end


BC=struct('abl',BCabl);

end

