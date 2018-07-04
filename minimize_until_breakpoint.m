function t = minimize_until_breakpoint(H, g, d, bp)

    dHd = d'*H*d;
    gd = g'*d;
    
    if gd > 0
        t = 0;
    elseif (gd <= 0 && dHd < 0) || (gd < 0 && dHd == 0)
        t = bp;
    else
        t = -gd/dHd;
        if t >= bp
            t = bp;
        end
    end

end