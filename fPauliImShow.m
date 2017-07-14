function z = fPauliImShow(varargin)

data = varargin{1};
n = 2;
if nargin ==2
    n = varargin{2};
end

z(:,:,3) = (data(:,:,1))./(mean2(data(:,:,1))*n);
z(:,:,1) = (data(:,:,2))./(mean2(data(:,:,2))*n);
z(:,:,2) = (data(:,:,3))./(mean2(data(:,:,3))*n);
