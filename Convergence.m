function [ out ] = OptimizeFinalError()
% Test script to check the convergence for different algorithms + some initial hyper parameter searching

algorithms = ["douglasrachfordprimal"];
parameters = ["gammal1", "rho", "t"];
images = ["testimages/cameraman.jpg"]; %"testimages/mcgill.jpg",, "testimages/manWithHat.tiff"
problems = ["l1", "l2"];
kernels = [struct("size", 5, 'sigma', 3), struct("size", 8, 'sigma', 5)]; % add noise


% gamma max, min, step size
% rho max, min, step size
% t max, min, step size
gammaMax = 0.2;
gammaMin = 0.05;
gammaStep = 0.05;

rhoMax = 0.5;
rhoMin = 0.1;
rhoStep = 0.1;

tMax = 4.0;
tMin = 1.0;
tStep = 1.0;


% Initialize output variable
out = struct('algorithm', [], 'gamma', [], 'rho', [], 't', [],  'error', [], 'kernels', [], 'images', [], 'problem', []);
% allocate enough room for all iterations, i.e out should have size #t diff
% * #gamma diff * etc...
out = repmat(struct('algorithm', '', 'gamma', [], 'rho', [], 't', [], 'error', [], 'kernels', [], 'images', '', 'problem', ''), ...
             numel(algorithms) * numel(parameters) * numel(images) * numel(problems) * numel(kernels) * ...
             numel(gammaMin:gammaStep:gammaMax) * numel(rhoMin:rhoStep:rhoMax) * numel(tMin:tStep:tMax), 1);
% Loop through algorithms, parameters, and images
totalnumberofiterations = numel(algorithms) * numel(parameters) * numel(images) * numel(problems) * numel(kernels)* ...
    numel(gammaMin:gammaStep:gammaMax) * numel(rhoMin:rhoStep:rhoMax) * numel(tMin:tStep:tMax);
iter=1;
for alg = algorithms
    for param = parameters
        for imgPath = images
            for problem = problems
                for kernelparams = kernels
                    for gamma = gammaMin : gammaStep : gammaMax
                        for rho = rhoMin : rhoStep : rhoMax
                            for t = tMin : tStep : tMax
                                img = PreprocessImage(imgPath);
                                if strcmp(imgPath, "testimages/mcgill.jpg")
                                    img = imresize(img, 0.05);
                                end
                                [numRows, numCols] = size(img);
                                kernel = fspecial("gaussian", kernelparams.size, kernelparams.sigma);
                                blurredimage = imfilter(img, kernel);
                                iterants = DefaultInitializeIterants(alg, blurredimage);
                                params.maxiter = 500;
                                if strcmp(problem, "l1")==1
                                    params.gammal1 = gamma;
                                else 
                                    params.gammal2 = gamma;
                                end
                                params.tprimaldr = t;
                                params.rhoprimaldr = rho;
                                params.blurredimage = blurredimage;
                                err = optsolve(problem, alg, iterants, kernel, blurredimage, params, img, true);
                                % out is the distance between unblurred and
                                % blurred
                                out(iter) = struct('algorithm', alg, 'gamma', gamma, 'rho', rho, 't', t, 'error', err, 'kernels', kernel, 'images', imgPath, 'problem', problem);
                                iter = iter + 1;
                                fprintf('Iter %4d | Time: %.2f completion \n', iter-1, ((iter-1)/totalnumberofiterations)*100);
                            end
                        end
                    end
                end
            end
        end
    end
end