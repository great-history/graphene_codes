function new_order = helper_get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, degeneracy, first_idxs, last_idxs, dims, eig_num, eps)
    % 分流：eigvec_H_m_K_now中简并度为1的本征态与eigvec_H_m_K_last中简并度为1的本征态有一一对应关系，即内积有且仅有一个为1，其余全为0
    % 而eigvec_H_m_K_now中多重简并态与eigvec_H_m_K_last中的多个本征态有交叠，即内积可能有多个不为0，但都小于1
    new_order = zeros(dims, 1);
    is_find = zeros(1,dims);
    for num = 1:eig_num
%         val_now = eigvals_unique(num);
        if degeneracy(num) == 1  % 简并度为1的态
            first = first_idxs(num);
            vec_now = eigvecs_now(:,first);
            val_now = eigvals_now(first);
            % 开始地毯式搜索
            target = 0;
            max_id = 0.0;
            for jj = 1:dims
                if is_find(jj) == 0 && abs(val_now - eigvals_last(jj)) < eps
                    inner_dot = abs(dot(vec_now, eigvecs_last(:,jj)));
                    if inner_dot > max_id  % 理想上应该是逼近于1的
                        max_id = inner_dot;
                        target = jj;
                    end
                end
            end
            new_order(target) = first;  % 注意不是new_order(first) = target !!!
            is_find(target) = 1;
            
        else  % 多重简并的态
            first = first_idxs(num);
            last = last_idxs(num);
            
            inner_dots = zeros(dims,1);
            for jj = 1:dims
                if is_find(jj) == 0
                    inner_dot = 0.0;
                    for kk = first:last
                        inner_dot = inner_dot + abs(dot(eigvecs_last(:,jj), eigvecs_now(:,kk)));
                    end
                    inner_dots(jj) = inner_dot;
                else
                    inner_dots(jj) = 0.0;
                end
            end
            
            [~, idxs] = sort(inner_dots);
            new_order(idxs(end-degeneracy(num)+1:end)) = first:last;
            is_find(idxs(end-degeneracy(num)+1:end)) = 1;
        end
    end
end