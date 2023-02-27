function [f, rec] = SplitBregAnalysis(g, A, AT, W, WT, opts)
%% SplitBregman Iterations for Solving
%   argmin 0.5*||Af-g||_2^2 + ||Wf||_1
%      f

    %% initialization
    if ~exist('opts','var') opts = struct; end
    if isfield(opts,'nIter') nIter = opts.nIter; else nIter = 30; end
    if isfield(opts,'nIterCG') nIterCG = opts.nIterCG; else nIterCG = 50; end
    if isfield(opts,'tolCG') tolCG = opts.tolCG; else tolCG = 1e-16; end
    if isfield(opts,'delta') delta = opts.delta; else delta = 1; end
    if isfield(opts,'mu') mu = opts.mu; else mu = 1; end
    if isfield(opts,'u0') u = opts.u0; else u = g; end
    if isfield(opts,'d0') d = opts.d0; else d = W(u); end
    if isfield(opts,'b0') b = opts.b0; else b = 0; end
    %% main loop
    CG_A = @(x) AT(A(x)) + mu*WT(W(x));
    for i = 1: nIter
        CG_B = AT(g) + mu*WT((d - b));
        u =  CG(CG_A, CG_B, u, tolCG, nIterCG);
        u(u>1) = 1; u(u<0) = 0;
        Wu = W(u);
        d = wthresh(Wu + b, 's', 1/mu);
        b = b + delta*(Wu - d);
        %%
        vecObj(i) = 0.5*norm(A(u)-g,'fro')^2 + sum(abs(Wu(:)),1);
        fprintf('iter=%d, obj_v=%f \n', i, vecObj(i));
    end
    f = u;
    rec.vecObj = vecObj;
end