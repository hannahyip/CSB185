% Read the data from an Excel file
filename = 'PreprocessedData.xlsx';
rpkm_data = readtable(filename);  % Read data into a table

genes = rpkm_data.Geneid;   % Convert cell array to string array
rpkm_data_values = rpkm_data{:, 2:end};  % Extract RPKM data (remove gene names column)
rpkm_data_values(:, end) = [];

% Perform k-means clustering on the RPKM data
numClusters = 2;
[clusterAssignments, ~] = kmeans(rpkm_data_values, numClusters, 'Replicates', 10);

% Perform PCA on RPKM data
[coeff, score, ~, ~, explained] = pca(rpkm_data_values);

% Plot first two principal components
figure;
gscatter(score(:,1), score(:,2), clusterAssignments, ['r', 'b'], 'o', 8);
xlabel(['Principal Component 1']);
ylabel(['Principal Component 2']);
title('Principal Component Analysis');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

% Pre-allocate a cell array for the table data
resultTable = cell(16564, 19); % 1 column for clusterAssignments + 18 columns for rpkm_data_values

% Populate the first column with clusterAssignments
resultTable(:, 1) = num2cell(clusterAssignments);

% Populate the remaining columns with rpkm_data_values
for i = 1:18
    resultTable(:, i+1) = num2cell(rpkm_data_values(:, i));
end

% Convert to a table
resultTable = cell2table(resultTable, 'VariableNames', ['ClusterAssignment', strcat('Gene_', string(1:18))]);
resultTableSorted = sortrows(resultTable, 'ClusterAssignment');
outputFile = 'output_file.xlsx';

% Write the result to a new Excel file
writetable(resultTableSorted, outputFile);

% Extract the gene expression data (excluding 'ClusterAssignment' column)
rpkm_data_values = resultTableSorted{:, 2:end};

% Calculate the median for rows 1 to 16142 (Cluster 1)
median_C1 = median(rpkm_data_values(1:16142, :), 1);

% Calculate the median for rows 16143 to 16564 (Cluster 2)
median_C2 = median(rpkm_data_values(16143:16564, :), 1);

% Display the median results
disp(median_C1);
disp(median_C2);
