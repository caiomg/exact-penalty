all_filenames = {
'20180403T193114_p1_db.mat';
'20180403T174453_p1_db.mat';
'20180403T175450_p1_db.mat';
'20180403T180837_p1_db.mat';
'20180403T181938_p1_db.mat';
'20180403T183008_p1_db.mat';
'20180403T183754_p1_db.mat';
'20180403T184828_p1_db.mat';
'20180403T190424_p1_db.mat';
'20180403T191600_p1_db.mat';
}



for k = 1:10
    name = all_filenames{k};
    load(name);
    df_results(k).mu = mu;
    df_results(k).results = results;
end
