function T = C2T(C)

% check argument
if nargin ~= 1
    error('Invalid argument.');
end

if ~isnumeric(C)
    error('Invalid argument.');
end

if ndims(C) ~= 3
    error('Invalid argument.');
end

[row, col, pag] = size(C);

if (row <= 0 || col <=0) || pag ~= 9
    error('Invalid argument.');
end

% convert to coherency matrix
T = zeros(row, col, 9);

T(:, :, 1) = (C(:, :, 1) + C(:, :, 3) + C(:, :, 5) * 2) / 2;
T(:, :, 2) = (C(:, :, 1) + C(:, :, 3) - C(:, :, 5) * 2) / 2;
T(:, :, 3) = C(:, :, 2);
T(:, :, 4) = (C(:, :, 1) - C(:, :, 3)) / 2;
T(:, :, 5) = (C(:, :, 4) + C(:, :, 6)) / sqrt(2);
T(:, :, 6) = (C(:, :, 4) - C(:, :, 6)) / sqrt(2);
T(:, :, 7) = -C(:, :, 8);
T(:, :, 8) = (C(:, :, 7) - C(:, :, 9)) / sqrt(2);
T(:, :, 9) = (C(:, :, 7) + C(:, :, 9)) / sqrt(2);
