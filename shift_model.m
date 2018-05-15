function model = shift_model(model, s)

    if isfield(model, 'c')
        model.c = model.c + model.g'*s + 0.5*(s'*model.H*s);
    end
    model.g = model.g + model.H*s;
    
end