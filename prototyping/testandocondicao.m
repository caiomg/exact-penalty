for k = 1:363
    if norm(hs2(k).x - solucao(k).x) > 0
        1;
    end
end