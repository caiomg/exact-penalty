function g_proj = l1_projected_gradient(p_grad, Q, R)

    [dim, r_cols] = size(R);
    if r_cols == 0
        g_proj = norm(p_grad);
    elseif r_cols < dim
        g_proj = norm(p_grad'*Q(:, r_cols+1:end));
    else
        g_proj = 0;
    end
    
end
