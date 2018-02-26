function r = rnd(alpha,beta,gamma,delta,varargin)
%STBL.RND   Alpha-stable random number generator
%   R = STBL.RND(ALPHA,BETA,GAMMA,DELTA) draws a sample from the Levy
%   alpha-stable distribution with characteristic exponent ALPHA, skewness
%   BETA, scale parameter GAMMA and location parameter DELTA. ALPHA, BETA,
%   GAMMA and DELTA must be size-compatible and fall in the following
%   ranges:
%
%         0 < ALPHA <= 2 (alpha=1 cannot be mixed with other cases)
%        -1 <= BETA <= 1
%         0 < GAMMA < inf
%      -inf < DELTA < inf
%
%   R = STBL.RND(ALPHA,BETA,GAMMA,DELTA,M,N,...) or
%   R = STBL.RND(ALPHA,BETA,GAMMA,DELTA,[M,N,...]) returns an M-by-N-by-...
%   array. Specified dimension must be size-compatible with the parameter 
%
%   References:
%   [1] J.M. Chambers, C.L. Mallows and B.W. Stuck (1976)
%       "A Method for Simulating Stable Random Variables"
%       JASA, Vol. 71, No. 354. pages 340-344
%
%   [2] Aleksander Weron and Rafal Weron (1995)
%       "Computer Simulation of Levy alpha-Stable Variables and Processes"
%       Lec. Notes in Physics, 457, pages 379-392

% Check parameters
narginchk(4,inf);
validateattributes(alpha,{'double','float'},{'positive','<=',2},mfilename,'ALPHA');
validateattributes(beta,{'double','float'},{'>=',-1,'<=',1},mfilename,'BETA');
validateattributes(gamma,{'double','float'},{'nonnegative','finite'},mfilename,'GAMMA');
validateattributes(delta,{'double','float'},{'finite'},mfilename,'DELTA');
if nargin==5 && ~isscalar(varargin{1})
   validateattributes(varargin{1},{'numeric'},{'row','integer','positive'},mfilename,'[M,N,...]')
elseif nargin>4
   cellfun(@(d)validateattributes(d,{'numeric'},{'scalar','integer','positive'},mfilename,'M,N,...'),...
      varargin);
end

if any(alpha==1) && ~all(alpha==1)
   error('stbl:rnd:invalidAlpha','ALPHA can either be all 1''s or all non-1''s.');
end

% Get the output size
sizeOut = getOutputSize(4,alpha,beta,gamma,delta,varargin{:});

%---Generate sample----

% See if parameters reduce to a special case, if so be quick, if not
% perform general algorithm

if all(alpha(:)==2)                          % Gaussian distribution
   r = sqrt(2) * randn(sizeOut);
   
elseif alpha(1)==1 && all(beta(:)==0)   % Cauchy distribution
   r = tan( pi/2 * (2*rand(sizeOut) - 1) );
   
elseif all(alpha(:) == .5) && all(abs(beta(:)) == 1) % Levy distribution (a.k.a. Pearson V)
   r = beta ./ randn(sizeOut).^2;
   
else
   V = pi/2 * (2*rand(sizeOut) - 1);
   W = -log(rand(sizeOut));
   if all(beta(:) == 0)        % Symmetric alpha-stable
      r = sin(alpha .* V) ./ ( cos(V).^(1./alpha) ) .* ...
         ( cos( V.*(1-alpha) ) ./ W ).^( (1-alpha)./alpha );
      
   elseif alpha(1) ~= 1   % General case, alpha not 1
      const = beta .* tan(pi/2*alpha);
      B = atan( const );
      S = (1 + const .* const).^(1./(2*alpha));
      r = S .* sin( alpha.*V + B ) ./ ( cos(V) ).^(1./alpha) .* ...
         ( cos( (1-alpha) .* V - B ) ./ W ).^((1-alpha)./alpha);
      
   else                        % General case, alpha = 1
      piover2 = pi/2;
      sclshftV =  piover2 + beta * V ;
      r = 1/piover2 * ( sclshftV .* tan(V) - ...
         beta * log( (piover2 * W .* cos(V) ) ./ sclshftV ) );
   end
end

% Scale and shift
r = gamma .* r + delta;
if alpha == 1 % Cauchy
   r(:) = r + (2/pi) * beta .* gamma .* log(gamma);
end
end

%====  function to find output size ======%
function commonSize = getOutputSize(nparams,varargin)

% set size vector based on 
if numel(varargin)==nparams
   commonSize = [1 1];
elseif iscell(varargin{nparams+1})
   commonSize = varargin{nparams+1}; % size vector given
elseif numel(varargin)==nparams+1
   commonSize = varargin{nparams+1}([1 1]);
else
   commonSize = [varargin{nparams+1:end}];
end

% expand as needed to match parameter array sizes
for argnum = 1:nparams
   tmp_size = size(varargin{argnum});
   tmp_ndims = ndims(varargin{argnum});
   commonSize(end+1:tmp_ndims) = 1; % expand if needed
   idx = commonSize(1:tmp_ndims)~=tmp_size;
   if any(commonSize(idx)~=1 & tmp_size(idx)~=1)
      error('stbl:rnd:inputSizeMismatch','Input argument sizes are not consistent.');
   end
   idx = tmp_size>1;
   commonSize(idx) = tmp_size(idx);
end
end
