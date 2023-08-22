% function H_hf_new = get_H_hf_by_diis(H_hf_diis_cell, error_vector_cell, num_error_vec, dims)
%     B_mat = zeros(num_error_vec + 1);
%     for ii = 1:num_error_vec
%         for jj = 1:num_error_vec
%             B_mat(ii, jj) = (error_vector_cell{ii})' * error_vector_cell{jj};
%         end
%     end
%     
%     for ii = 1:num_error_vec
%         B_mat(num_error_vec + 1, ii) = -1;
%         B_mat(ii, num_error_vec + 1) = -1;
%     end
%     
%     C_vec = zeros(num_error_vec + 1, 1);
%     C_vec(end) = -1;
%     X_coeff_list = linsolve(B_mat, C_vec);
%     
%     H_hf_new = zeros(dims);
%     for ii = 1:num_error_vec
%         H_hf_new = H_hf_new + X_coeff_list(ii) * H_hf_diis_cell{ii};
%     end
%     
% end


function H_hf_new = get_H_hf_by_diis(H_hf_diis_cell, error_vector_cell, num_error_vec, dims)
    % error_vector_cell中装的是
    B_mat = zeros(num_error_vec + 1);
    for ii = 1:num_error_vec
        for jj = 1:num_error_vec
            B_mat(ii, jj) = trace(error_vector_cell{ii} * (error_vector_cell{jj})');
        end
    end
    
    for ii = 1:num_error_vec
        B_mat(num_error_vec + 1, ii) = -1;
        B_mat(ii, num_error_vec + 1) = -1;
    end
    
    C_vec = zeros(num_error_vec + 1, 1);
    C_vec(end) = -1;
    X_coeff_list = linsolve(B_mat, C_vec);
    
    H_hf_new = zeros(dims);
    for ii = 1:num_error_vec
        H_hf_new = H_hf_new + X_coeff_list(ii) * H_hf_diis_cell{ii};
    end
    
end