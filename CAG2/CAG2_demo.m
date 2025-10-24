clear;
clc;
addpath(genpath('./'));
datadir='./datasets/';
dataname={'proteinFold'};
numdata = length(dataname);
numname = {''};

% 设置运行次数
num_runs = 20;

for idata = 1:length(dataname)
    ResBest = zeros(1, 8);
    ResStd = zeros(1, 8);
    
    % Initialize results table
    allResults = []; 
    colNames = {'m', 'lambda', 'beta', ...
                'ACC', 'NMI', 'Purity', 'F1', 'Precision', 'Recall', ...
                'ARI', 'Fowlkes', 'runtime'};
    
    for dataIndex = 1:1
        datafile = [datadir, dataname{idata}, numname{dataIndex}, '.mat'];
        load(datafile);
        if exist('Y', 'var') && exist('X', 'var')
            truelabel = {Y};
            gt = Y;
            data = X;
            n = size(gt, 1);
            numview = length(data);
         
            % Transpose each view to d_v x n format
            for v = 1:numview
                data{v} = data{v}'; 
            end
            
            % Create complete data index matrix
            ind = true(n, numview);
            X1 = data;
        end
        k = length(unique(gt));
        numview = length(X1);
        n = size(X1{1}, 2); 
        %% proteinFold
        TempLambda1 = [5*k];  
        TempLambda2 = [0.001];
        TempLambda3 = [1e-6];
        
        
        
        % Initialize storage for parameter grid results
        ACC = zeros(length(TempLambda1), length(TempLambda2), length(TempLambda3));
        NMI = zeros(length(TempLambda1), length(TempLambda2), length(TempLambda3));
        Purity = zeros(length(TempLambda1), length(TempLambda2), length(TempLambda3));
        Runtime = zeros(length(TempLambda1), length(TempLambda2), length(TempLambda3));
        idx = 1;
        
        % Main
        for LambdaIndex1 = 1:length(TempLambda1)
            m_val = TempLambda1(LambdaIndex1);
            for LambdaIndex2 = 1:length(TempLambda2)
                lambda_val = TempLambda2(LambdaIndex2);
                for LambdaIndex3 = 1:length(TempLambda3)
                    beta_val = TempLambda3(LambdaIndex3);
                    % Display current parameter configuration
                    fprintf('\n%s - m=%d, lambda=%.4f, beta=%.4f\n', dataname{idata}, m_val, lambda_val, beta_val);
                    
                    % ===== MULTIPLE RUNS FOR EACH PARAMETER COMBINATION =====
           
                    
                    % Initialize storage for multiple runs
                    all_ACC = zeros(num_runs, 1);
                    all_NMI = zeros(num_runs, 1);
                    all_Purity = zeros(num_runs, 1);
                    all_F1 = zeros(num_runs, 1);
                    all_Precision = zeros(num_runs, 1);
                    all_Recall = zeros(num_runs, 1);
                    all_ARI = zeros(num_runs, 1);
                    all_Fowlkes = zeros(num_runs, 1);
                    all_runtime = zeros(num_runs, 1);
                    all_obj_val = zeros(num_runs, 1);
                    
                    % Multiple runs loop
                    for run = 1:num_runs
                        
                        % Set random seed for reproducibility
                        rng(run);
                        
                        tic;
                        % Run optimization
                        [Zall, P, A, Z_client, obj_val] = CAG2_opt(X1, ind, m_val, lambda_val, beta_val, 50);
                        
                        % Compute SVD
                        [UU,~,~] = svd(Zall', 'econ');
                        UU = UU(:, 1:k);
                        F = UU ./ sqrt(sum(UU.^2, 2));
                        
                        % Handle NaN values
                        nan_idx = any(isnan(F), 2);
                        if any(nan_idx)
                            F = F(~nan_idx, :);
                            gt_clean = gt(~nan_idx);
                        else
                            gt_clean = gt;
                        end
                        
                        % Run k-means
                        pY = kmeans(F, k, 'maxiter', 1000, 'replicates', 10, 'emptyaction', 'singleton');
                        
                        % Compute metrics
                        run_results = Clustering8Measure(gt_clean, pY);
                        runtime = toc;
                        
                        % Store results
                        all_ACC(run) = run_results(1);
                        all_NMI(run) = run_results(2);
                        all_Purity(run) = run_results(3);
                        all_F1(run) = run_results(4);
                        all_Precision(run) = run_results(5);
                        all_Recall(run) = run_results(6);
                        all_ARI(run) = run_results(7);
                        all_Fowlkes(run) = run_results(8);
                        all_runtime(run) = runtime;
                        all_obj_val(run) = obj_val;
                    end
                    
                    % Compute mean and standard deviation
                    tempResBest = [mean(all_ACC), mean(all_NMI), mean(all_Purity), mean(all_F1), ...
                                  mean(all_Precision), mean(all_Recall), mean(all_ARI), mean(all_Fowlkes)];
                    tempResStd = [std(all_ACC), std(all_NMI), std(all_Purity), std(all_F1), ...
                                 std(all_Precision), std(all_Recall), std(all_ARI), std(all_Fowlkes)];
                    runtime_mean = mean(all_runtime);
                    runtime_std = std(all_runtime);
                    
                    % Print statistics
                    fprintf('    ACC: %.4f ± %.4f\n', tempResBest(1), tempResStd(1));
                    
                    
                    % Update grid results
                    ACC(LambdaIndex1, LambdaIndex2, LambdaIndex3) = tempResBest(1);
                    NMI(LambdaIndex1, LambdaIndex2, LambdaIndex3) = tempResBest(2);
                    Purity(LambdaIndex1, LambdaIndex2, LambdaIndex3) = tempResBest(3);
                    Runtime(LambdaIndex1, LambdaIndex2, LambdaIndex3) = runtime_mean;
                    
                    % Build row for Excel output
                    newRow = [m_val, lambda_val, beta_val, ...
                              tempResBest(1), tempResBest(2), tempResBest(3), ...
                              tempResBest(4), tempResBest(5), tempResBest(6), ...
                              tempResBest(7), tempResBest(8), runtime_mean];
                    allResults = [allResults; newRow];
                    
                    % Update global best results (based on mean ACC)
                    if tempResBest(1) > ResBest(1, 1)
                        newZ = Zall;
                        newF = F;
                        ResBest = tempResBest;
                        ResStd = tempResStd;
                    end
                    idx = idx + 1;
                end
            end
        end
        
        % Save aggregated results (removed)
        aRuntime = mean(Runtime(:));
        PResBest = ResBest;
        PResStd = ResStd;
    end
    
    % Export comprehensive Excel file (removed)
    
    % Save final aggregated results (removed)
    

end