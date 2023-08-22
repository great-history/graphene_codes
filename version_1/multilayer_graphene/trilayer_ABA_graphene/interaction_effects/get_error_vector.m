function error_vector = get_error_vector(H_hf, density_mat, dims, type)
    % 首先计算对易
    error_vector = H_hf * density_mat - density_mat * H_hf;
    % 然后把error vector拉直
    if type == "vector"
        error_vector = reshape(error_vector, [dims * dims, 1]);
    end
end