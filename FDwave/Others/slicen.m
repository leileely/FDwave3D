function hout = slicen(varargin)
% A new version of the slice plot.
%   SLICE(X,Y,Z,V,Sx,Sy,Sz) draws slices along the x,y,z directions at
%   the points in the vectors Sx,Sy,Sz. The arrays X,Y,Z define the
%   coordinates for V and must be monotonic and 3-D plaid (as if
%   produced by MESHGRID).  The color at each point will be determined
%   by 3-D interpolation into the volume V.  V must be an M-by-N-by-P
%   volume array. 
%
%   SLICE(X,Y,Z,V,XI,YI,ZI) draws slices through the volume V along the
%   surface defined by the arrays XI,YI,ZI.
%
%   SLICE(V,Sx,Sy,Sz) or SLICE(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P. 
%
%   SLICE(...,'method') specifies the interpolation method to use.
%   'method' can be 'linear', 'cubic', or 'nearest'.  'linear' is the
%   default (see INTERP3).
%
%   SLICE(AX,...) plots into AX instead of GCA. 
%
%   H = SLICE(...) returns a vector of handles to SURFACE objects.
%
%   The axes CLim property is set to span the finite values of V.
%
%   Example: To visualize the function x*exp(-x^2-y^2-z^2) over the
%   range -2 < x < 2, -2 < y < 2, -2 < z < 2, 
%
%      [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%      v = x .* exp(-x.^2 - y.^2 - z.^2);
%      slice(x,y,z,v,[-1.2 .8 2],2,[-2 -.2])
%
%   See also MESHGRID, INTERP3.

%   J.N. Little 1-23-92
%   Copyright 1984-2010 The MathWorks, Inc. 

%%   06/03/2019
%%   a small bug for homogeneous media in line 140 is fixed based on slice.m in Matlab 2018a.

[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(4,8,nargs,'struct'));

cax = newplot(cax);
hold_state = ishold(cax);
axis_slice = 0; % Slice only along axes?
method = 'linear';

nin = nargs;
if ischar(args{nin}), % slice(...,'method')
  method = args{nin};
  nin = nin - 1;
end

if nin == 4 || nin == 5 % slice(v,xi,yi,zi,<nx>)
    v = args{1}; args{1} = [];

    if nin == 5,
      warning(message('MATLAB:slice:ExtraShortArgument'))
      if ndims(v)==2, % Reshape the volume.
        v = reshape(v,size(v,1)/args{5},args{5},size(v,2));
      end
    end

    if ndims(v)~=3, error(message('MATLAB:slice:InvalidVDims')); end
    [ny,nx,nz] = size(v);
    if min(size(args{2}))<=1 && min(size(args{3}))<=1 && ...
       min(size(args{4}))<=1,
        axis_slice = 1;
    end
    [sx,sy,sz] = deal(args{2:4});
    [x,y,z] = meshgrid(1:nx,1:ny,1:nz);
elseif nin == 6,
    error(message('MATLAB:slice:WrongNumberOfInputs'));
elseif nin == 7 || nin == 8, % slice(x,y,z,v,xi,yi,zi,<nx>)
    v = args{4}; args{4} = [];

    if nin == 8,
      warning(message('MATLAB:slice:ExtraArgument'))
      if ndims(v)==2, % Reshape the volume
        v = reshape(v,size(v,1)/args{8},args{8},size(v,2));
      end
    end

    if ndims(v)~=3, error(message('MATLAB:slice:InvalidVDims')); end
    if min(size(args{1}))==1 || min(size(args{2}))==1 || ...
       min(size(args{3}))==1,
        [x,y,z] = meshgrid(args{1:3});
    else
        [x,y,z] = deal(args{1:3});
    end
    [ny,nx,nz] = size(v);

    if min(size(args{5}))<=1 && min(size(args{6}))<=1 && ...
       min(size(args{7}))<=1,
        axis_slice = 1;
    end
    [sx,sy,sz] = deal(args{5:7});
end

if axis_slice,
    h = [];
    [xi,yi,zi] = meshgrid(sx,y(:,1,1),z(1,1,:));
    vi = interp3(x,y,z,v,xi,yi,zi,method);
    for i = 1:length(sx)
        h = [h; surface( ...
            reshape(xi(:,i,:),[ny nz]),reshape(yi(:,i,:),[ny nz]), ...
            reshape(zi(:,i,:),[ny nz]),reshape(vi(:,i,:),[ny nz]), ...
            'parent',cax)];
    end

    [xi,yi,zi] = meshgrid(x(1,:,1),sy,z(1,1,:));
    vi = interp3(x,y,z,v,xi,yi,zi,method);
    for i = 1:length(sy)
        h = [h; surface( ...
            reshape(xi(i,:,:),[nx nz]),reshape(yi(i,:,:),[nx nz]),...
            reshape(zi(i,:,:),[nx nz]),reshape(vi(i,:,:),[nx nz]), ...
            'parent',cax)];
    end

    [xi,yi,zi] = meshgrid(x(1,:,1),y(:,1,1),sz);
    vi = interp3(x,y,z,v,xi,yi,zi,method);
    for i = 1:length(sz)
        h = [h; surface(xi(:,:,i),yi(:,:,i),zi(:,:,i),vi(:,:,i), ...
            'parent',cax)];
    end
else
    vi = interp3(x,y,z,v,sx,sy,sz,method);
    h = surf(sx,sy,sz,vi,'parent',cax);
end

if nargout > 0
    hout = h;
end
if ~hold_state
    view(cax,3), grid(cax,'on')
end

% Use ISFINITE to make sure no NaNs or Infs get passed to CAXIS
u=v(isfinite(v)); u = u(:);
% raw
%caxis(cax,[min(u) max(u)]);
% new
clim = [min(u) max(u)];
if(clim(1) ~= clim(2))
    caxis(cax,clim);
end
% Signal to the world that we have created a new plot:
plotdoneevent(cax,h);


