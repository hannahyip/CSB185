% Read the data from an Excel file
filename = 'PreprocessedData.xlsx';
rpkm_data = readtable(filename);  % Read data into a table

genes = rpkm_data.Geneid;   % Convert cell array to string array
rpkm_data_values = rpkm_data{:, 2:end};  % Extract RPKM data (remove gene names column)
rpkm_data_values(:, end) = [];

% Fit the Gaussian mixture model with 2 clusters
numClusters = 2;
gm = fitgmdist(rpkm_data_values, numClusters, 'Replicates', 5, 'Options', statset('MaxIter', 500));

% Assign each data point to a cluster
clusterLabels = cluster(gm, rpkm_data_values);

% Save clustering results
outputFile = 'clustered_data.xlsx';

outputTable = array2table([rpkm_data_values, clusterLabels], 'VariableNames', [...
    strcat('Feature', string(1:size(rpkm_data_values,2))), 'Cluster']);
outputTableSorted = sortrows(outputTable, 'Cluster');
writetable(outputTableSorted, 'clustered_data.xlsx');

disp('Clustering complete. Results saved to clustered_data.xlsx.');

% Extract the gene expression data (excluding 'ClusterAssignment' column)
rpkm_data_values = outputTableSorted{:, [1:18, 20:end]};

% Calculate the median for rows 1 to 9849 (Cluster 1)
median_C1 = median(rpkm_data_values(1:9849, :), 1);

% Calculate the median for rows 9850 to 16564 (Cluster 2)
median_C2 = median(rpkm_data_values(9850:16564, :), 1);

% Display the median results
disp(median_C1);
disp(median_C2);

% Perform PCA on the RPKM data (excluding the cluster column)
[coeff, score, ~, ~, explained] = pca(rpkm_data_values);

% Select the first two principal components for visualization
pc1 = score(:, 1);
pc2 = score(:, 2);

% Scatter plot of PCA results with cluster colors
figure;
gscatter(pc1, pc2, outputTableSorted.Cluster, ['r', 'b'], 'o');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('Principal Component Analysis');
legend('Cluster 1', 'Cluster 2');
grid on;
