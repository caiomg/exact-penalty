if exist('terminate_cutest_problem', 'file') ~= 2
  addpath('my_cutest_functions');
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



% Parameters
mu = 10;
epsilon = 2;
delta = 1e-6;
Lambda = 0.1;

final_filenames = {};
for iter = 1:10
    clear results;

    fprintf(log_fd, ['mu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                     'Lambda = % 8g\n\n'], mu, epsilon, delta, Lambda);
    fprintf(1, ['mu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                     'Lambda = % 8g\n\n'], mu, epsilon, delta, Lambda);


    n_problems = length(all_problems);
    for k = 1:n_problems
        problem_name = all_problems{k};
        prob = setup_cutest_problem(problem_name, '../my_problems/');

        % Objective
        f_obj = @(x) get_cutest_objective(x);
        counter = evaluation_counter(f_obj);
        f = @(x) counter.evaluate(x);

        % Constraints
        n_constraints = get_cutest_total_number_of_constraints();
        all_con = cell(n_constraints, 1);
        for n = 1:n_constraints
            gk = @(x) evaluate_my_cutest_constraint(x, n, 1);
            all_con{n} = gk;
        end

        % Initial point
        x0 = prob.x;

        nlcon = @(x) constraints(all_con, {}, x, 1);
        fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                       'SpecifyObjectiveGradient', true);
        x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
        fcount_fmincon = counter.get_count();
        fx_fmincon = f(x_fmincon);
        nphi_fmincon = norm(max(0, nlcon(x_fmincon)));
        %%


        counter.reset_count();
        counter.set_max_count(5000);

        solved = true;
        try
            [x, hs2] = l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda);
        catch error
            results(k, 1).except = error;
            solved = false;
        end
        if solved
            fcount = counter.get_count();
            fx = f(x);
        nphi = norm(max(0, nlcon(x)));
            error_obj = fx_fmincon - fx;
            error_x = norm(x_fmincon - x);
        else
            x = [];
            hs2 = [];
            fx = [];
        nphi = [];
            error_obj = [];
            error_x = [];
            fcount = [];
        end

        results(k, 1).name = problem_name;
        results(k, 1).x = x;
        results(k, 1).fx = fx;
        results(k, 1).history = hs2;
        results(k, 1).fcount = fcount;
        results(k, 1).fcount_fmincon = fcount_fmincon;
        results(k, 1).error_obj = error_obj;
        results(k, 1).error_x = error_x;
        results(k, 1).nphi = nphi;
        results(k, 1).nphi_fmincon = nphi_fmincon;

        print_results(results(k, 1), log_fd);
        print_results(results(k, 1));
        filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
        save(filename, 'results', 'mu', 'epsilon', 'delta', 'Lambda', 'final_filenames');

        terminate_cutest_problem();
    end
    final_filenames{iter, 1} = filename;
    mu = mu*2;
end
