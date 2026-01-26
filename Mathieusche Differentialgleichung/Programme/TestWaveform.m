close all; clear; clc;
load('all_data_matrix_w03.mat', 'all_data_matrix' )

load('results_by_branch_w03.mat','results_by_branch' )

index0 = all_data_matrix(:,3) == 0; 
all_data_matrix0 = all_data_matrix(index0, :);
figure
plot(all_data_matrix0 (: ,1),all_data_matrix0(:,2));
index1 = all_data_matrix(:,3) == -1;
all_data_matrix1 = all_data_matrix(index1, :);
hold on; plot(all_data_matrix1 (: ,1),all_data_matrix1(:,2));