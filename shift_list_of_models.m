function models_shifted = shift_list_of_models(models, s)
% SHIFT_LIST_OF_MODELS - 

    for k = numel(models):-1:1
        models_shifted(k) = shift_model(models(k), s);
    end
end
