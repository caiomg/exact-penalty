statuswarning = warning('query', 'MATLAB:declareGlobalBeforeUse');
warning('off', 'MATLAB:declareGlobalBeforeUse');
global xg
global fmodel
global myconstraints

if exist('gfx', 'var')
    fmodel.H = Hfx;
    fmodel.g = gfx;
    xg = x;
    myconstraints = phi;
end
warning(statuswarning.state, 'MATLAB:declareGlobalBeforeUse');