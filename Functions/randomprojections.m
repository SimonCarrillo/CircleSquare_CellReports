function [X_new] = randomprojections(X,n_out)

n_samples = size(X,1);
n_features = size(X,2);
n_iterations = 100;

X_tmp = nan(n_samples,n_out,n_iterations);

for iter = 1:n_iterations
    
    projected_matrix_tmp = rand(n_features,n_out);
    projected_matrix = normalize(projected_matrix_tmp,'zscore');

    X_tmp(:,:,iter) = X*projected_matrix;
    
end

X_new = mean(X_tmp,3);

end