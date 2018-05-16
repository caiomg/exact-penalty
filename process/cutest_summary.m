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

n_problems = length(all_problems);

tol_c = 1e-3;
tol_g = 1e-3;

sif_directory = '/home/caio/local/cutest/sif';
sol_str = '*LO SOLTN';
for n = 1:n_problems
    problem = all_problems{n};
    filename = fullfile(sif_directory, [problem, '.SIF']);
    sif_contents = fileread(filename);
    ind_sol = regexp(sif_contents, sol_str, 'end') + 1;
    solution = sscanf(sif_contents(ind_sol:end), '%f');
    result_sif(n, 1).solution = solution;
end            
            
package = 'cobyla';
r_directory = fullfile('..', 'my_problems');


result_cutest(n_problems, 1).f_count = [];
result_cutest(n_problems, 1).fval = [];
for n = 1:n_problems
    problem = all_problems{n};
    this_directory = fullfile(r_directory, problem);
    summary_filename = [problem '_result_' package];
    summary_filename = fullfile(this_directory, summary_filename);
    [fc, fv, x] = read_cutest_summary(summary_filename);
    result_cutest(n, 1).name = problem;
    result_cutest(n, 1).f_count = fc;
    result_cutest(n, 1).fval = fv;
    result_cutest(n, 1).x = x;
    result_cutest(n, 1).sol_fval = result_sif(n, 1).solution;
end

for n = 1:n_problems
    problem = all_problems{n};
    terminate_cutest_problem();
    prob = setup_cutest_problem(problem, '../my_problems/');

    % Objective
    f_obj = @(x) get_cutest_objective(x);

    % Constraints
    n_constraints = get_cutest_total_number_of_constraints();
    all_con = cell(n_constraints, 1);
    for k = 1:n_constraints
        gk = @(x) evaluate_my_cutest_constraint(x, k, 1);
        all_con{k} = gk;
    end
    nlcon = @(x) constraints(all_con, {}, x, 1);

    if ~isempty(result_cutest(n).sol_fval)
        result_cutest(n).error = result_cutest(n).fval - result_cutest(n).sol_fval;
    end
    if ~isempty(result_cutest(n, 1).x)
        result_cutest(n).real_c = norm(max(0, nlcon(result_cutest(n, 1).x)));
        result_cutest(n).real_f = f_obj(result_cutest(n, 1).x);
        [is_kkt, lgrad] = check_kkt(f_obj, all_con, result_cutest(n, 1).x, tol_c, tol_g);
        result_cutest(n).kkt = is_kkt;
        result_cutest(n).lgrad = lgrad;
    end
    terminate_cutest_problem();
end

results_cobyla = result_cutest

