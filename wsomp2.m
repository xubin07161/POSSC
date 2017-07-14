function gamma = wsomp2(varargin)
%WSOMP2 Error-constrained Weighted Simultaneous Orthogonal Matching Pursuit.
%  GAMMA = WSOMP2(D,X,W,LISTGROUP,EPSILON) solves the optimization problem
%
%       min  |GAMMA_i|_0  s.t. |W.*(X_i - D*GAMMA_i)|_F <= n*n_i*EPSILON^2
%      gamma
%
%  for each of the groups in X, using Weighted Orthogonal Matching Pursuit.
%  Here, D is a dictionary with normalized columns, X is a matrix
%  containing column signals, which may be grouped in various sizes,
%	 and W is the weight assigned to each entry of X.
%  EPSILON is the error target for each entry,
%  The output GAMMA is a matrix containing
%  the sparse representations as its columns. 
%
%
%  GAMMA = WSOMP2(D,X,W,DtXWW,DDtWW,LISTGROUP,EPSILON) is the fastest implementation of OMP2,
%  but also requires the most memory. Here, DtXWW stores the weighted projections
%  D'*(X.*W.*W), and DDtWW is a matrix containing the squared norms of the
%  weighted dictionary, (D.*D)'*(W.*W). 
%  Note that in general, the call
%
%    GAMMA = WSOMP2(D,X,1./X_std,D'*(X./X_std.^2),(D.*D)'*(1./X_std.^2), LISTGROUP, EPSILON);
%
%  will be faster than the call
%
%    GAMMA = WSOMP2(D,X,1./X_std,LISTGROUP,EPSILON);
%
%  due to optimized matrix multiplications in Matlab. However, when the
%  entire matrix D'*(X.*W.*W) cannot be stored in memory, one of the other two
%  versions can be used. Both compute D'*(X.*W.*W) for just one group at a time,
%  and thus require much less memory.
%
%  GAMMA = WSOMP2(...,[],EPSILON) suggests that all group have the same
%  group size specified in maxgroupsize, in this case, both 'groupnumber'
%  and 'maxgroupsize' should be specified.
%
%  GAMMA = WSOMP2(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies additional
%  parameters for WSOMP2. Available parameters are:
%	   'groupnumber'  - Specified the number of groups
%	
%    'maxgroupsize' - Specified the maximun size of group, when list group
%                     is not specified, this implies uniform group size
%
%    'gammamode' - Specifies the representation mode for GAMMA. Can be
%                  either 'full' or 'sparse', corresponding to a full or
%                  sparse matrix, respectively. By default, GAMMA is
%                  returned as a sparse matrix.
%    'maxatoms' -  Limits the number of atoms in the representation of each
%                  signal. If specified, the number of atoms in each
%                  representation does not exceed this number, even if the
%                  error target is not met. Specifying maxatoms<0 implies
%                  no limit (default).
%    'messages'  - Specifies whether progress messages should be displayed.
%                  When positive, this is the number of seconds between
%                  status prints. When negative, indicates that no messages
%                  should be displayed (this is the default).
%    'profile'   - Can be either 'on' or 'off'. When 'on', profiling
%                  information is displayed at the end of the funciton
%                  execution.
%
%
%  Summary of WSOMP2 versions:
%
%    version																		|   speed     |   memory
%  ------------------------------------------------------------------------
%   WSOMP2(D,X,W,DtXWW,DDtWW,LISTGROUP,EPSILON) |  fast				|  large
%   WSOMP2(D,X,W,LISTGROUP,EPSILON)							|  slow				|  small
%  ------------------------------------------------------------------------
%
%
%  References:
%  [1] M. Elad, R. Rubinstein, and M. Zibulevsky, "Efficient Implementation
%      of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit",
%      Technical Report - CS, Technion, April 2008.
%
%  See also SOMP.

%  Original author:
%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009

%  Modified by:
%  Bin Zuo
%  Electronic Engineering Department
%  Tsinghua University, Beijing, 100084, P.R.China
%  zuob2009@126.com

%  January 2014


% default options

S = 0;
maxl = 0;
sparse_gamma = 1;
msgdelta = -1;
maxatoms = -1;
profile = 0;


% determine number of parameters

paramnum = 1;
while (paramnum<=nargin && ~ischar(varargin{paramnum}))
  paramnum = paramnum+1;
end
paramnum = paramnum-1;


% parse options

for i = paramnum+1:2:length(varargin)
  paramname = varargin{i};
  paramval = varargin{i+1};

  switch lower(paramname)
		case 'groupnumber'
			S = paramval;
			
		case 'maxgroupsize'
			maxl = paramval;

    case 'gammamode'
      if (strcmpi(paramval,'sparse'))
        sparse_gamma = 1;
      elseif (strcmpi(paramval,'full'))
        sparse_gamma = 0;
      else
        error('Invalid GAMMA mode');
      end
      
    case 'maxatoms'
      maxatoms = paramval;

    case 'messages'
      msgdelta = paramval;

    case 'profile'
      if (strcmpi(paramval,'on'))
        profile = 1;
      elseif (strcmpi(paramval,'off'))
        profile = 0;
      else
        error('Invalid profile mode');
      end

    otherwise
      error(['Unknown option: ' paramname]);
  end
  
end


% determine call type

if (paramnum==5)	%  D,X,W,LISTGROUP,EPSILON
     
	D = varargin{1};
	X = varargin{2};
	W = varargin{3};
	DtXWW = [];
	DDtWW = [];
	listgroup = varargin{4};
	epsilon = varargin{5};  
elseif (paramnum == 7) % D,X,W,DtXWW,DDtWW,LISTGROUP,EPSILON
	D = varargin{1};
	X = varargin{2};
	W = varargin{3};
	DtXWW = varargin{4};
	DDtWW = varargin{5};
	listgroup = varargin{6};
	epsilon = varargin{7}; 
else
	error('Invalid number of parameters');
end


% verify dictionary normalization
% 
% if (isempty(G))
%   atomnorms = sum(D.*D);
% else
%   atomnorms = diag(G);
% end
% if (any(abs(atomnorms-1) > 1e-2))
%   error('Dictionary columns must be normalized to unit length');
% end

% verify legality of listgroup
if ~isempty(listgroup)
	if size(listgroup,1) > 1 && size(listgroup,2) > 1
		error('Listgroup must be a vector');
	end
	if listgroup(1) == 1 % if listgroup begins with 1
		listgroup = listgroup - 1;
	end
	if listgroup(end) < size(X,2)
		SS = length(listgroup);
		mmaxl = max(diff([listgroup(:);size(X,2)]));
		if S && maxl && (S ~= SS || mmaxl ~= maxl)
			warning('Group number and group size is inconsistent with listgroup, the former two are ignored');
		end
		S = SS;
		maxl = mmaxl;
	else
		error('List_group have elements too large.');
	end
elseif S && maxl
	if S*maxl ~= size(X,2)
		error('Group number and group size are incompatible with X.');
	end
else
	error('Either group number and group size, or list_group, must be specified.');
end

% WSomp

gamma = wsomp2mex(D,X,W,DtXWW,0,DDtWW,listgroup,epsilon,S,maxl,sparse_gamma,msgdelta,maxatoms,profile);

end