classdef GridClass < handle
    %GRIDCLASS defines grid for velocity models
    % Authors are temporarily removed for review.
    % Publication date: 13th of July, 2019
    properties
        x0
        dx
        nx
        y0
        dy
        ny
        z0
        dz
        nz
        t0
        dt
        nt
        mx
        my
        mz
        mt
        xx
        yy
        zz
        tt
    end
    
    methods
        function valid = checkGrid(obj)
            if obj.nx ~=0
                valid = 1; 
            else 
                valid = 0; 
            end
        end
        %
        function gridInfo(obj)
            disp('Information about grid:'); 
            disp(['x0=' num2str(obj.x0) ', dx=' num2str(obj.dx) ', Nx=' num2str(obj.nx) '.']);
            disp(['y0=' num2str(obj.y0) ', dy=' num2str(obj.dy) ', Ny=' num2str(obj.ny) '.']);
            disp(['z0=' num2str(obj.z0) ', dz=' num2str(obj.dz) ', Nz=' num2str(obj.nz) '.']);
            disp(['t0=' num2str(obj.t0) ', dt=' num2str(obj.dt) ', Nt=' num2str(obj.nt) '.']);   
        end
        %
        function setGrid(obj)
            obj.mx = (obj.nx-1)*obj.dx + obj.x0; 
            obj.my = (obj.ny-1)*obj.dy + obj.y0; 
            obj.mz = (obj.nz-1)*obj.dz + obj.z0; 
            obj.mt = (obj.nt-1)*obj.dt + obj.t0; 

            obj.xx = obj.x0:obj.dx:obj.mx; 
            obj.yy = obj.y0:obj.dy:obj.my; 
            obj.zz = obj.z0:obj.dz:obj.mz; 
            obj.tt = obj.t0:obj.dt:obj.mt; 
        end
        %
        function G = oldGrid(obj)
            G(1) = obj.x0;
            G(2) = obj.y0; 
            G(3) = obj.z0; 
            G(4) = obj.nx; 
            G(5) = obj.ny; 
            G(6) = obj.nz; 
            G(7) = obj.dx; 
            G(8) = obj.dy; 
            G(9) = obj.dz; 
            G(10) = obj.t0; 
            G(11) = obj.nt; 
            G(12) = obj.dt; 
        end
        %
        function newGrid(obj, Gold)
            obj.x0 = Gold(1);
            obj.y0 = Gold(2); 
            obj.z0 = Gold(3); 
            obj.nx = Gold(4); 
            obj.ny = Gold(5); 
            obj.nz = Gold(6); 
            obj.dx = Gold(7); 
            obj.dy = Gold(8); 
            obj.dz = Gold(9); 
            obj.t0 = Gold(10); 
            obj.nt = Gold(11); 
            obj.dt = Gold(12); 
        end
    end
end


%% for C++

% GridClass G;



% class GridClass{
%     public:
%         int    nx;
%         int    ny;
%         int    nz;
%         int    nt;
%         double x0;
%         double dx;
%         double y0;
%         double dy;
%         double z0;
%         double dz;
%         double t0;
%         double dt;
%         double mx;
%         double my;
%         double mz;
%         double mt;
%         double xx;
%         double yy;
%         double zz;
%         double tt;
% };


% // Get G - parameters of the velocity grid
%     if (!mxIsClass(prhs[2], "GridClass")) mexErrMsgTxt("Input (5th arg.) must be an object of class GridClass");
%     if (mxGetProperty(prhs[2],0,"nx")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'nx' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"nx"));
% 	G.nx = int(*tmp);
%     if (mxGetProperty(prhs[2],0,"ny")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'ny' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"ny"));
% 	G.ny = int(*tmp);
%     if (mxGetProperty(prhs[2],0,"nz")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'nz' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"nz"));
% 	G.nz = int(*tmp);
%     if (mxGetProperty(prhs[2],0,"nt")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'nt' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"nt"));
% 	G.nt = int(*tmp);
%     if (mxGetProperty(prhs[2],0,"x0")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'x0' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"x0"));
% 	G.x0 = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"y0")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'y0' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"y0"));
% 	G.y0 = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"z0")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'z0' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"z0"));
% 	G.z0 = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"t0")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 't0' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"t0"));
% 	G.t0 = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"dx")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'dx' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"dx"));
% 	G.dx = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"dy")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'dy' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"dy"));
% 	G.dy = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"dz")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'dz' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"dz"));
% 	G.dz = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"dt")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'dt' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"dt"));
% 	G.dt = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"mx")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'mx' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"mx"));
% 	G.mx = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"my")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'my' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"my"));
% 	G.my = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"mz")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'mz' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"mz"));
% 	G.mz = double(*tmp);
%     if (mxGetProperty(prhs[2],0,"mt")==NULL)   mexErrMsgTxt("Input (5th arg): Required Property 'mt' is missing.");
% 	tmp = (double *)mxGetPr(mxGetProperty(prhs[2],0,"mt"));
% 	G.mt = double(*tmp);