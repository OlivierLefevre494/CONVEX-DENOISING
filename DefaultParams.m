function params = DefaultParams()
params.maxiter = 500;

params.gammal1 = 0.01;
params.gammal2 = 0.05;

params.tprimaldr = 0.3;
params.tprimaldualdr = 16;

params.rhoprimaldualdr = 0.8;
params.rhoprimaldr = 0.8;

params.tadmm   = 16;
params.rhoadmm = 1.5;

params.tchamb = 0.35;
params.schamb = 0.35; %Highly dependant on A thus harder to set one specific case
end

% Note that these are picked from the hyperparam table and picked somewhat
% from experience as no single hyperparameter can hit all problem types by
% default