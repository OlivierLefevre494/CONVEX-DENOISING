function iterants = DefaultInitializeIterants(algorithm)
if strcmp(algorithm, "douglasrachfordprimal")==1
    iterants.z1 = zeros(size(blurredimage)) + 0.5;
    iterants.z2 = zeros([size(blurredimage), 3]) + 0.5;
elseif strcmp(algorithm, "douglasrachfordprimaldual")
    iterants.p = zeros(size(blurredimage)) + 0.5;
    iterants.q = zeros([size(blurredimage), 3]) + 0.5;
elseif strcmp(algorithm, "admm")==1
    % Add initialization here
elseif strcmp(algorithm, "chambollepock")==1
    % Add initialization here
end
end