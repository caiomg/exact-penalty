function [s, change] = truncated_cg_step(H, g, radius)

dim = size(g, 1);
gk = g;
s = zeros(dim, 1);
pk = -gk;

for k = 1:dim
    php = pk'*H*pk;
    if php <= 0
        break
    end
    a = (gk'*gk)/php;
    sn = s + a*pk;
    if norm(sn) > radius || a == 0
       % Stop iterating
       break
    end
    gn = gk + a*(H*pk);
    b = (gn'*gn)/(gk'*gk);
    pk = -gn + b*pk;
    gk = gn;
    s = sn;
end

if norm(pk, inf) == 0
    1;
end
    
if norm(pk, inf) ~= 0 && (php <= 0 || norm(sn) > radius)
    % We have to truncate step at tr border
    a = roots([pk'*pk, 2*s'*pk, s'*s - radius^2]);
    a = max(a);
    if isempty(a) || ~isreal(a) || a < 0
       error();
    end
    s = s + a*pk;
end

change = g'*s + 0.5*(s'*H*s);

end