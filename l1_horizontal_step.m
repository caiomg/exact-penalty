function h = l1_horizontal_step(B, g, N, x0, radius, bl, bu)

    if isempty(bl)
        bl = -inf(size(x0));
    end
    if isempty(bu)
        bu = inf(size(x0));
    end

    % Evaluating matrices on reduced space
    Br = N'*B*N;
    gr = N'*g;

    % One shot of More-Sorensen algorithm
    s0 = ms_step(Br, gr, radius);
    h = N*s0;
    % One shot of projected-gradient step
%     h = pg_path_bounds(B, g, x0, h0, bl, bu, radius);

end

