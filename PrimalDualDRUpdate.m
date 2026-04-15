function [ out ] = PrimalDualDRUpdate(iterants, problem, blurredimage, params, ApplyA, invertMatrix, ApplyATrans, last)
if last
    out = NormalizeImage(iterants.p);
else

    iterants.x = NormalizeImage(iterants.p); % This is the prox of f
    if strcmp(problem, "l2")==1 % scaled prox for l1 and isonorm
        iterants.z = cat(3, params.tprimaldualdr .* TwoNormProx(iterants.q(:,:,1) ./params.tprimaldualdr, blurredimage, 1/params.tprimaldualdr), params.tprimaldualdr .* IsoProx(iterants.q(:,:,2)./(params.tprimaldualdr), iterants.q(:,:,3)./(params.tprimaldualdr), params.gammal2/(params.tprimaldualdr)));
    elseif strcmp(problem, "l1")==1 % scaled prox for l2 and isonorm
        iterants.z = cat(3, params.tprimaldualdr .* OneNormProx(iterants.q(:,:,1) ./params.tprimaldualdr, blurredimage, 1/params.tprimaldualdr), params.tprimaldualdr .* IsoProx(iterants.q(:,:,2)./(params.tprimaldualdr), iterants.q(:,:,3)./(params.tprimaldualdr), params.gammal1/(params.tprimaldualdr)));
    end

%     if strcmp(problem, "l2")==1 % scaled prox for l1 and isonorm
%         iterants.z = cat(3, params.tprimaldualdr .* TwoNormProx(iterants.q(:,:,1), blurredimage./params.tprimaldualdr, 1/params.tprimaldualdr), params.tprimaldualdr .* IsoProx(iterants.q(:,:,2)./(params.tprimaldualdr), iterants.q(:,:,3)./(params.tprimaldualdr), params.gammal2/(params.tprimaldualdr)));
%     elseif strcmp(problem, "l1")==1 % scaled prox for l2 and isonorm
%         iterants.z = cat(3, params.tprimaldualdr .* OneNormProx(iterants.q(:,:,1), blurredimage./params.tprimaldualdr, 1/params.tprimaldualdr), params.tprimaldualdr .* IsoProx(iterants.q(:,:,2)./(params.tprimaldualdr), iterants.q(:,:,3)./(params.tprimaldualdr), params.gammal1/(params.tprimaldualdr)));
%     end

    % prox for conjugate of l1 and isonorm (we use Moreau identity)
    iterants.z = iterants.q - iterants.z;
    
    iterants.w = invertMatrix(2.0 .* iterants.x - iterants.p - params.tprimaldualdr .* ApplyATrans(2.0 .* iterants.z - iterants.q));
    iterants.v = params.tprimaldualdr .* ApplyA(iterants.w) + 2.0 .* iterants.z - iterants.q;
    iterants.p = iterants.p + params.rhoprimaldualdr .* (iterants.w - iterants.x);
    iterants.q = iterants.q + params.rhoprimaldualdr .* (iterants.v - iterants.z);
    
    
    out = iterants;
end
end