function  [s, fs, ind_eactive] = cauchy_step(m0, radius, N, mu, con_models, Ii0, d0, w0, ind_eactive, epsilon)
%CAUCHY_STEP Summary of this function goes here
%   Detailed explanation goes here

if isempty(d0)
    d0 = -m0.g;
end
d = d0;

if isempty(w0)
    w = zeros(size(d0));
else
    w = w0;
end

Ii = Ii0;

n_constraints = length(con_models);
gamma = [];
IM = [];

if isempty(N) || isempty(d) || norm(d) == 0
    s = w;
    fs = 0;
    return
end

% Calculating TR constraint
bpoint = roots([d'*d, 2*w'*d, w'*w - radius^2]);
gamma_tr = min(bpoint(bpoint > 0));
if isempty(gamma_tr) || ~isreal(gamma_tr)
    error;
end



for n = 1:n_constraints
    if ~isempty(find(Ii == n, 1))
        if con_models(n).c ~= 0
            Hc = con_models(n).H;
            gc = con_models(n).g;
            cc = con_models(n).c;
            bpoint = roots([0.5*(d'*(N'*Hc*N)*d), w'*(N'*Hc*N)*d + gc'*N*d, 0.5*w'*(N'*Hc*N)*w + gc'*N*w + cc]);
            if size(bpoint, 1) >= 1 && bpoint(1) > 0 && isreal(bpoint(1)) && bpoint(1) < gamma_tr
                gamma(end+1, 1) = bpoint(1);
                IM(end+1, 1) = n;
            end
            if size(bpoint, 1) >= 2 && bpoint(2) > 0 && isreal(bpoint(2)) && bpoint(2) < gamma_tr
                gamma(end+1, 1) = bpoint(2);
                IM(end+1, 1) = n;
            end
        end
    end
end
[gamma, ind_temp] = sort(gamma);
IM = IM(ind_temp);


B = m0.B;
g = m0.g;

tl = 0;
tu = Inf;
alpha = 0;
total_descent = 0;
constraint_changed = 0;

while ~isempty(IM)
    tu = gamma(1);
    n = IM(1);
    gamma = gamma(2:end);
    IM = IM(2:end);
    
    if (B'*w + g)'*d > 0
        % Ascent
        m1 = m0;
        m1.B = B;
        m1.g = g;
        w1 = w + alpha*d;
        if constraint_changed > 0
            nc = constraint_changed;
        elseif constraint_changed < 0
            nc = -constraint_changed;
        else
            break;
        end
        d1 = d - N'*(con_models(nc).H*N*w1 + con_models(nc).g)*...
                 (d'*N'*(con_models(nc).H*N*w1 + con_models(nc).g))/...
                 ((con_models(nc).H*N*w1 + con_models(nc).g)'*...
                 (con_models(nc).H*N*w1 + con_models(nc).g));
        A1 = [];
        for e = 1:n_constraints
            if (isempty(find(Ii == e, 1)) && ~isempty(find(Ii0 == e, 1))) || e == nc
                A1 = [A1, (con_models(nc).H*N*w1 + con_models(nc).g)];
            end
        end
        A1 = N'*A1;
        d1 = correct_direction(d1, A1);
        [s, segment_descent] = cauchy_step(m1, radius, N, mu, con_models, Ii, d1, w1, ind_eactive, epsilon);
        pstep = w1 - s;
        if segment_descent > mu*(0.5*(pstep'*N'*con_models(nc).H*N*pstep) + (con_models(nc).H'*N*w1 + con_models(nc).g)'*N*pstep)
            % Accept
            total_descent = total_descent + segment_descent;
            if constraint_changed < 0
                total_descent = total_descent - mu*(0.5*(pstep'*N'*con_models(nc).H*N*pstep) + (con_models(nc).H'*N*w1 + con_models(nc).g)'*N*pstep);
            end
        else
            % Rollback
            s = w1;
        end
        fs = total_descent;
        break;
    elseif d'*B*d < 0
        segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*tu^2 + (B'*w + g)'*d*tu);
        total_descent = total_descent + segment_descent;
        alpha = tu;
%         segment_descent = (0.5*(d'*B*d)*alpha^2 + g'*d*alpha) - (0.5*(d'*B*d)*tu^2 + g'*d*tu);
%         total_descent = total_descent + segment_descent;
%         alpha = tu;
    else
        % Minimizing the quadratic model in direction d
        tau = -((d'*B*d)\((B'*w + g)'*d));
        if tau < alpha
            % does this occur?
            s = w + alpha*d;
            fs = total_descent;
            break;
        elseif tau < tu
            segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*tau^2 + (B'*w + g)'*d*tau);
            total_descent = total_descent + segment_descent;
            alpha = tau;
            s = w + alpha*d;
            fs = total_descent;
            % change the active set?
            break
        else
            segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*tu^2 + (B'*w + g)'*d*tu);
            total_descent = total_descent + segment_descent;
            alpha = tu;
        end
    end
    if sign(alpha*((d'*N')*con_models(n).H*(N*d)) + con_models(n).g'*N*d) > 0
        B = B + mu*N'*con_models(n).H*N;
        g = g + mu*N'*con_models(n).g;
        constraint_changed = n;
    else
        B = B - mu*N'*con_models(n).H*N;
        g = g - mu*N'*con_models(n).g;
        constraint_changed = -n;
    end
    if isempty(find(IM == n, 1))
        Ii = Ii(Ii ~= n);
    end
end
if ~exist('s', 'var')
    
    if (B'*w + g)'*d > 0
        % Ascent
        s = w + alpha*d;
        fs = total_descent;
    elseif d'*B*d < 0
        segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*gamma_tr^2 + (B'*w + g)'*d*gamma_tr);
        alpha = gamma_tr;
        s = w + alpha*d;
        total_descent = total_descent + segment_descent;
        fs = total_descent;
    else
        % Minimizing the quadratic model in direction d
        tau = -((w'*B)*d + g'*d)/(d'*B*d);
        if tau < alpha
            % does this occur?
            s = w + alpha*d;
            fs = total_descent;
        elseif tau < gamma_tr
            segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*tau^2 + (B'*w + g)'*d*tau);
            total_descent = total_descent + segment_descent;
            alpha = tau;
            s = w + alpha*d;
            fs = total_descent;
        else
            segment_descent = (0.5*(d'*B*d)*alpha^2 + (B'*w + g)'*d*alpha) - (0.5*(d'*B*d)*gamma_tr^2 + (B'*w + g)'*d*gamma_tr);
            total_descent = total_descent + segment_descent;
            alpha = gamma_tr;
            s = w + alpha*d;
            fs = total_descent;
        end
    end
end
new_active = new_tr_eactive_constraints(con_models, Ii0, N*s, epsilon);
ind_eactive = [ind_eactive;
               new_active];


end

