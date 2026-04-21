function [ out ] = PrimalDRUpdate(iterants, problem, blurredimage, params, ApplyA, invertMatrix, ApplyATrans, last)
if last
    out = NormalizeImage(iterants.z1);
else

    iterants.x = NormalizeImage(iterants.z1);
    if strcmp(problem, "l2")==1
        iterants.y = cat(3, TwoNormProx(iterants.z2(:,:,1), blurredimage, params.tprimaldr), IsoProx(iterants.z2(:,:,2), iterants.z2(:,:,3), params.gammal2*params.tprimaldr));
    elseif strcmp(problem, "l1")==1
        iterants.y = cat(3, OneNormProx(iterants.z2(:,:,1), blurredimage, params.tprimaldr), IsoProx(iterants.z2(:,:,2), iterants.z2(:,:,3), params.gammal1*params.tprimaldr));
    end
    
    iterants.b = 2.*iterants.y-iterants.z2;
    iterants.u = invertMatrix(double(2.*iterants.x-iterants.z1)+ApplyATrans(iterants.b));
    iterants.v = ApplyA(iterants.u);
    iterants.z1 = double(iterants.z1) + (params.rhoprimaldr).*(double(iterants.u)-double(iterants.x));
    iterants.z2 = double(iterants.z2) + (params.rhoprimaldr).*(double(iterants.v)-double(iterants.y));
    
    out = iterants;
    
end
iterants.output = iterants.x;
end