filename = 'ClusteringMedians.xlsx'; % Clustering medians saved in excel file

% Read IFN-beta and IFN-lambda data
ifn_beta = xlsread(filename, 'B:L'); % Columns B-L
ifn_lambda = xlsread(filename, 'N:X'); % Columns N-X

% Extract Cluster 1 (Row 2) and Cluster 2 (Row 3)
cluster1_beta = ifn_beta(1, :)'; % Cluster 1 for IFN-beta (column format)
cluster2_beta = ifn_beta(2, :)'; % Cluster 2 for IFN-beta (column format)
cluster1_lambda = ifn_lambda(1, :)'; % Cluster 1 for IFN-lambda (column format)
cluster2_lambda = ifn_lambda(2, :)'; % Cluster 2 for IFN-lambda (column format)

% Save to .mat files
save('cluster1_beta.mat', 'cluster1_beta');
save('cluster2_beta.mat', 'cluster2_beta');
save('cluster1_lambda.mat', 'cluster1_lambda');
save('cluster2_lambda.mat', 'cluster2_lambda');