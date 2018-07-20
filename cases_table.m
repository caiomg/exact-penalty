if exist('terminate_cutest_problem', 'file') ~= 2
  addpath('my_cutest_functions');
end
if exist('check_kkt', 'file') ~= 2
  addpath('process');
end
if exist('evaluate_polynomial', 'file') ~= 2
  addpath('tr_modeling');
  addpath('tr_modeling/polynomials');
end



terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

logdir = fullfile('..', 'logs', datestr(now, 30));
if ~exist(logdir, 'dir')
   if mkdir(logdir) ~= 1
       logdir = '.';
   end
end
log_filename = fullfile(logdir, sprintf('%s_p1_db.log', datestr(now, 30)));
log_fd = fopen(log_filename, 'w');


    all_problems = {'BURKEHAN' 'HIMMELP2' 'HIMMELP3' 'HIMMELP4' 'HIMMELP5' ...
                    'HIMMELP6' 'HS10' 'HS11' 'HS12' 'HS13' 'HS15' 'HS16' ...
                    'HS17' 'HS18' 'HS19' 'HS20' 'HS21' 'HS22' 'HS23' ...
                    'HS24' 'HS57' 'HS59' 'HS88' 'HUBFIT' 'LSQFIT' ...
                    'SIMPLLPA' 'SIMPLLPB' 'SNAKE' 'TWOBARS' 'ZECEVIC2' ...
                    'ZECEVIC3' 'ZECEVIC4' 'CB2' 'CB3' 'CHACONN1' ...
                    'CHACONN2' 'CONGIGMZ' 'DEMYMALO' 'GIGOMEZ1' 'GIGOMEZ2' ...
                    'GIGOMEZ3' 'HS29' 'HS30' 'HS31' 'HS33' ...
                    'HS34' 'HS35' 'HS35I' 'HS35MOD' 'HS36' 'HS37' 'HS64' ...
                    'HS65' 'HS66' 'HS67' 'HS89' 'KIWCRESC' 'LOOTSMA' ...
                    'MADSEN' 'MAKELA1' 'MAKELA2' 'MIFFLIN1' 'MIFFLIN2' ...
                    'MINMAXRB' 'POLAK1' 'POLAK4' 'POLAK5' 'SPIRAL' ...
                    'STANCMIN' 'WOMFLET' 'ZY2' 'BIGGSC4'  ...
                    'HATFLDH' 'HS43' 'HS44' 'HS44NEW' 'HS70' 'HS72' ...
                    'HS76' 'HS76I' 'HS90' 'CANTILVR' 'EXPFITA' 'HS268' ...
                    'HS83' 'HS84' 'HS85' 'HS86' 'HS91' 'MINMAXBD' 'POLAK6' ...
                    'ROSENMMX' 'S268' 'CRESC4' 'HALDMADS' 'HS92' 'HS93' ...
                    'HS95' 'HS96' 'HS97' 'HS98' 'MATRIX2' 'PENTAGON' ...
                    'SYNTHES1' 'DIPIGRI' 'HS100' 'HS100MOD' 'HS101' ...
                    'HS102' 'HS103' 'HS21MOD' 'S365' 'S365MOD' 'AVGASA' ...
                    'AVGASB' 'HS104' 'HS105' 'HS106' 'EQC' 'HS108' ...
                    'MISTAKE' 'QC' 'QCNEW' 'HS113' 'POLAK2' 'POLAK3' ...
                    'HAIFAS' 'HS116' 'HS117' 'HS118' 'DEMBO7' 'MAKELA3' ...
                    'MAKELA4' 'OPTPRLOC' };
                

l1_options = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 0.5, 'radius_max', 1e3, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'diagonal hessian', ...
                        'pivot_threshold', 0.1, 'poised_radius_factor', 2, ...
                        'pivot_imp', 1.1)

% Parameters
mu = 500;
epsilon = 0.85;
delta = 1e-6;
Lambda = 0.075;

list_of_problems

all_mu = [10, 50, 100, 1000, 10000, 100000, 1000000]

clear tries

all_epsilon = [1]
all_lambda = [0.075]

tries(1).epsilon = 0.85;
tries(1).Lambda = 0.075; %best

warning('off', 'cmg:bad_fvalue');

final_filenames = {};
all_results = {};
good_results = {};
all_solved = [];
iter = 0;
n_problems = length(selected_problems);
solved_problems = false(n_problems, 1);
for mu_i = 1:length(all_mu)
    for la_i = 1:length(all_lambda)
        for ep_i = 1:length(all_epsilon)
            iter = iter + 1;
            epsilon = all_epsilon(ep_i);
            Lambda = all_lambda(la_i);
            mu = all_mu(mu_i);

            clear results;

            fprintf(log_fd, ['\nmu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                             'Lambda = % 8g\n'], mu, epsilon, delta, Lambda);
            fprintf(1, ['\nmu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                             'Lambda = % 8g\n'], mu, epsilon, delta, Lambda);


            for k = 1:n_problems
                if solved_problems(k)
                    continue
                end
                problem_name = selected_problems(k).name;
                prob = setup_cutest_problem(problem_name, '../my_problems/');

                % Objective
                f_obj = @(x) get_cutest_objective(x);
                counter = evaluation_counter(f_obj);
                f = @(x) counter.evaluate(x);

                % Constraints
                n_constraints = get_cutest_total_number_of_constraints();

                bl = [];
                bu = [];
                % Bound constraints
%                 bl = prob.bl;
                % bu = prob.bu;
                lower_bounds = prob.bl > -1e19;
                upper_bounds = prob.bu < 1e19;
                cutest_lower_bounds = prob.bl > -1e19;
                cutest_upper_bounds = prob.bu < 1e19;
                lower_bounds = bl > -1e19;
                upper_bounds = bu < 1e19;

                % Other constraints
                all_con = cell(prob.m, 1);

                for p = 1:(n_constraints - sum(cutest_lower_bounds) - ...
                           sum(cutest_upper_bounds))
                    gk = @(x) evaluate_my_cutest_constraint(x, p, 1);
                    all_con{p} = gk;
                end
                % Lower bounds
                c_ind = p;
                for p = 1:prob.n
                    if cutest_lower_bounds(p)
                        c_ind = c_ind + 1;
                        if isempty(lower_bounds) || ~lower_bounds(p)
                            gk = @(x) evaluate_my_cutest_constraint(x, c_ind, 1);
                            all_con{end+1} = gk;
                        end
                    end
                end
                % Upper bounds
                c_ind = length(all_con);
                for p = 1:prob.n
                    if cutest_upper_bounds(p)
                        c_ind = c_ind + 1;
                        if isempty(upper_bounds) || ~upper_bounds(p)
                            gk = @(x) evaluate_my_cutest_constraint(x, c_ind, 1);
                            all_con{end+1} = gk;
                        end
                    end
                end
                if length (all_con) ~= n_constraints - sum(lower_bounds) ...
                        - sum(upper_bounds)
                    error();
                end


                % Initial point
                x0 = prob.x;

                nlcon = @(x) constraints(all_con, {}, x, 1);
%                 fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
%                                                'SpecifyObjectiveGradient', true);
%                 [x_fmincon, fx_fmincon, exitflag, output, mult_fmincon] = fmincon(f, x0,[],[],[],[],bl,bu, nlcon, fmincon_options);
%                 fprintf(1, '| %8s |', problem_name);
%                 fprintf(1, '  % +9.3g |\n', sum(abs([mult_fmincon.lower; mult_fmincon.upper; mult_fmincon.ineqnonlin])));

%                 fcount_fmincon = counter.get_count();
%                 fx_fmincon = f(x_fmincon);
%                 nphi_fmincon = norm(max(0, nlcon(x_fmincon)));
                %%

                counter.reset_count();
                counter.set_max_count(30000);

                try
                    p_seed = rng('default');
                    [x, hs2] = l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, bl, bu, l1_options);
                    solved = true;
                catch thiserror
                    results(k, 1).except = thiserror;
                    solved = false;
                end
                rng(p_seed);
                if solved
                    fcount = counter.get_count();
                    fx = f(x);
                    nphi = norm(max(0, nlcon(x)));
                    error_obj = selected_problems(k).solution - fx;
                    [kkt, lgrad] = check_kkt(f, all_con, x, bl, bu, 1e-6, 1e-4);
                else
                    x = [];
                    hs2 = [];
                    fx = [];
                    nphi = [];
                    error_obj = [];
                    error_x = [];
                    fcount = nan;
                    kkt = false;
                    lgrad = nan;
                end

                results(k, 1).name = problem_name;
                results(k, 1).x = x;
                results(k, 1).fx = fx;
                results(k, 1).history = []; % hs2
                results(k, 1).fcount = fcount;
                results(k, 1).error_obj = error_obj;
                results(k, 1).nphi = nphi;
                results(k, 1).mu = mu;
                results(k, 1).epsilon = epsilon;
                results(k, 1).lambda = Lambda;
                results(k, 1).kkt = kkt;
                results(k, 1).lgrad = lgrad;

                print_results(results(k, 1), log_fd);
                print_results(results(k, 1));

                all_results{iter} = results;
                if kkt || ...
                        (~isempty(nphi) && nphi < 1e-6 && abs(error_obj) < 5e-6)
                    solved_problems(k) = true;
                    good_results{sum(solved_problems)} = results(k, 1);
                    all_solved(end+1) = k;
                end
                terminate_cutest_problem();
            end
                filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
                save(filename, 'all_results', 'mu', 'epsilon', ...
                     'delta', 'Lambda', 'final_filenames', 'good_results', 'all_solved');
        end
    end
end
