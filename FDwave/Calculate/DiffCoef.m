function C = DiffCoef(o,s)
% The function is used to calculate differential coefficient
% o: the half precision order of difference
% s: grid scheme
% 'r' indicates regular-grid scheme;'s' indicates stagged-grid scheme
% Examples:
%       C = DCoef(o) return 2*o order coefficient of regular-grid scheme
%       C = DCoef(o,'r') the same as above
%       C = DCoef(o,'s') return 2*o order coefficient of stagged-grid scheme
C = zeros(1,o)';A = zeros(o);
if nargin < 1
    error('Please enter input arguments!');
elseif nargin == 1
    B = [1/2,zeros(1,o-1)]';
    for i=1:o
        for j=1:o
            A(i,j)=j^(2*i-1);
        end
    end
elseif nargin>2
    error('Two many input arguments!');
else
    if s=='r'
        B = [1/2,zeros(1,o-1)]';
        for i=1:o
            for j=1:o
                A(i,j)=j^(2*i-1);
            end
        end
    elseif s=='s'
        B = [1,zeros(1,o-1)]';
        for i=1:o
            for j=1:o
                A(i,j)=(2*j-1)^(2*i-1);
            end
        end
    else
        error('''s'' can only be ''r'' or ''s''!');
    end
end
C=A\B;