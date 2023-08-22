% 单层石墨烯克服Dirac point的奇异性
gamma0 = 3.1;
dims = 2;
Nx = 200;
Ny = 200;
err = 0.2;
akxs = linspace(-0.1, 0.1, Nx);
akys = linspace(-0.1, 0.1, Ny);

eig_enes_K = zeros(Nx, Ny, 2);
eig_enes_Kp = zeros(Nx, Ny, 2);

for i = 1:Nx
    akx = akxs(i);
    aky = akys(1);
    [H_m_K, H_m_Kp] = construct_monolayer_effective_Ham(0, 0, gamma0, akx, aky);

    [eigvec_H_m_K_last, D] = eig(H_m_K);
    eigval_H_m_K_diag_last = diag(D);

    [eigvec_H_m_Kp_last, D] = eig(H_m_Kp);
    eigval_H_m_Kp_diag_last = diag(D);

    eig_enes_K(i,1,:)=eigval_H_m_K_diag_last;
    eig_enes_Kp(i,1,:)=eigval_H_m_Kp_diag_last;

    for j = 2:Ny
        aky = akys(j);
        [H_m_K, H_m_Kp] = construct_monolayer_effective_Ham(0, 0, gamma0, akx, aky);

        [eigvec_H_m_K_now, D] = eig(H_m_K);
        eigval_H_m_K_diag_now = diag(D);

        [eigvec_H_m_Kp_now, D] = eig(H_m_Kp);
        eigval_H_m_Kp_diag_now = diag(D);

        [eigvec_H_m_K_now, eigval_H_m_K_diag_now, eig_num_K, eigvec_H_m_Kp_now, eigval_H_m_Kp_diag_now, eig_num_Kp] = ...
            re_order_monolayer(eigvec_H_m_K_last, eigval_H_m_K_diag_last, eigvec_H_m_K_now, eigval_H_m_K_diag_now, ...
                               eigvec_H_m_Kp_last, eigval_H_m_Kp_diag_last, eigvec_H_m_Kp_now, eigval_H_m_Kp_diag_now, ...
                               err);

        eig_enes_K(i, j, :)=eigval_H_m_K_diag_now;
        eig_enes_Kp(i, j, :)=eigval_H_m_Kp_diag_now;

        % 如果全是非简并，那么我们需要将这一次的结果，即eigvec_H_m_K_now/eigval_H_m_K_diag_now和eigvec_H_m_Kp_now/eigval_H_m_Kp_diag_now用到下一次
        % 如果存在简并的情形，那么我们需要将上一次的结果，即eigvec_H_m_K_last/eigval_H_m_K_diag_last和eigvec_H_m_Kp_last/eigval_H_m_Kp_diag_last用到下一次
        if (eig_num_K == dims)
           eigvec_H_m_K_last = eigvec_H_m_K_now;
           eigval_H_m_K_diag_last = eigval_H_m_K_diag_now;
        end

        if (eig_num_Kp == dims)
           eigvec_H_m_Kp_last = eigvec_H_m_Kp_now;
           eigval_H_m_Kp_diag_last = eigval_H_m_Kp_diag_now;
        end
    end
end

[aakkxs, aakkys]=meshgrid(akxs, akys);

figure % 1
mesh(aakkxs,aakkys,eig_enes_K(:,:,1))
hold on
mesh(aakkxs,aakkys,eig_enes_K(:,:,2))
colorbar

figure % 2
subplot(1,2,1);
for i = 1:Nx
    plot(akxs,eig_enes_K(i,:,1),'LineWidth', 2)
    hold on;
    plot(akxs,eig_enes_K(i,:,2),'LineWidth', 2)
    hold on;
end

grid on;
xlabel('k_x \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-1 -1/2 0 1/2 1]);
xticklabels({'-1' '-1/2' '0' '1/2' '1'})
title('effective model E_k versus k_x; for k_y = 0 at K valley')
legend('Valence band', 'Conduction band');
hold on

subplot(1,2,2);
plot(akxs,eig_enes_Kp(101,:,1),'LineWidth', 2)
grid on;
hold on;
plot(akxs,eig_enes_Kp(101,:,2),'LineWidth', 2)
xlabel('k_x \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-1 -1/2 0 1/2 1]);
xticklabels({'-1' '-1/2' '0' '1/2' '1'})
title('effective model E_k versus k_x; for k_y = 0 at Kp valley')
legend('Valence band', 'Conduction band');



function [eigvec_H_m_K_now, eigval_H_m_K_diag_now, eig_num_K, eigvec_H_m_Kp_now, eigval_H_m_Kp_diag_now, eig_num_Kp] = ...
            re_order_monolayer(eigvec_H_m_K_last, eigval_H_m_K_diag_last, eigvec_H_m_K_now, eigval_H_m_K_diag_now, ...
                               eigvec_H_m_Kp_last, eigval_H_m_Kp_diag_last, eigvec_H_m_Kp_now, eigval_H_m_Kp_diag_now, ...
                               err)
    % 要保证eigval_H_m_K_last和eigval_H_m_Kp_last都是非简并的才行！！！
    dims = 2;
    % 首先检查eigval_H_m_K_last和eigval_H_m_Kp_last是否都是非简并的
    [eigvals_K_unique, degeneracy_K, eig_num_K, first_idxs_K, last_idxs_K] = get_degeneracy(eigval_H_m_K_diag_now, dims);
    [eigvals_Kp_unique, degeneracy_Kp, eig_num_Kp, first_idxs_Kp, last_idxs_Kp] = get_degeneracy(eigval_H_m_Kp_diag_now, dims);
    
    % 分流：eigvec_H_m_K_now中简并度为1的本征态与eigvec_H_m_K_last中简并度为1的本征态有一一对应关系，即内积有且仅有一个为1，其余全为0
    % 而eigvec_H_m_K_now中多重简并态与eigvec_H_m_K_last中的多个本征态有交叠，即内积可能有多个不为0，但都小于1
    % is_find = zeros(dims,1);
    
    new_order_K = get_new_order(eigvec_H_m_K_last, eigval_H_m_K_diag_last, eigvec_H_m_K_now, eigvals_K_unique, degeneracy_K, first_idxs_K, last_idxs_K, dims, eig_num_K, err);
    new_order_Kp = get_new_order(eigvec_H_m_Kp_last, eigval_H_m_Kp_diag_last, eigvec_H_m_Kp_now, eigvals_Kp_unique, degeneracy_Kp, first_idxs_Kp, last_idxs_Kp, dims, eig_num_Kp, err);
    
     % 交换位置
    eigvec_H_m_K_now(:,:) = eigvec_H_m_K_now(:,new_order_K);
    eigval_H_m_K_diag_now(:) = eigval_H_m_K_diag_now(new_order_K);
    
    eigvec_H_m_Kp_now(:,:) = eigvec_H_m_Kp_now(:,new_order_Kp);
    eigval_H_m_Kp_diag_now(:) = eigval_H_m_Kp_diag_now(new_order_Kp);
end

function [eigvals_unique, degeneracy, eig_num, first_idxs, last_idxs] = get_degeneracy(eigvals, dims)
    % 假定eigval_H_m_K_diag_last是从小到大排列的，所以返回的简并度可以很方便地转换到指标，
    % 比如[1,1,2,3,1] → [1,2,4,7,8]就是本征值的第一位置
    degeneracy = zeros(dims,1);
    eigvals_unique = zeros(dims,1);
    
    eig_val_now = eigvals(1);
    eig_num = 1;
    eigvals_unique(eig_num) = eig_val_now;
    count = 1;
    
    for i = 2:dims
        eig_val_last = eig_val_now;
        eig_val_now = eigvals(i);
        if eig_val_now == eig_val_last
            count = count + 1;
        else
            degeneracy(eig_num) = count;
            eig_num = eig_num + 1;
            eigvals_unique(eig_num) = eig_val_now;
            count = 1;
        end
    end
    degeneracy(eig_num) = count;
    
    % pull out elements
    eigvals_unique = eigvals_unique(1:eig_num);
    degeneracy = degeneracy(1:eig_num);
    
    first_idxs = zeros(eig_num,1);
    last_idxs = zeros(eig_num,1);
    
    last = 0;
    for num = 1:eig_num
        first = last + 1;
        last = last + degeneracy(num); 
        first_idxs(num) = first;
        last_idxs(num) = last;
    end
end

function new_order = get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_unique, degeneracy, first_idxs, last_idxs, dims, eig_num, err)
    % 分流：eigvec_H_m_K_now中简并度为1的本征态与eigvec_H_m_K_last中简并度为1的本征态有一一对应关系，即内积有且仅有一个为1，其余全为0
    % 而eigvec_H_m_K_now中多重简并态与eigvec_H_m_K_last中的多个本征态有交叠，即内积可能有多个不为0，但都小于1
    new_order = zeros(dims,1);
    for num = 1:eig_num
%         val_now = eigvals_unique(num);
        if degeneracy(num) == 1  % 简并度为1的态
            first = first_idxs(num);
            vec_now = eigvecs_now(:,first);
            
            % 开始地毯式搜索
            target = 0;
            max_id = 0.0;
            for jj = 1:dims
                inner_dot = abs(dot(vec_now, eigvecs_last(:,jj)));
                if inner_dot > max_id  % 理想上应该是逼近于1的
                    max_id = inner_dot;
                    target = jj;
                end
            end
            new_order(first) = target;
            
        else  % 多重简并的态
            first = first_idxs(num);
            last = last_idxs(num);
            
            inner_dots = zeros(dims,1);
            for jj = 1:dims
                inner_dot = 0.0;
                for kk = first:last
                    inner_dot = inner_dot + abs(dot(eigvecs_last(:,jj), eigvecs_now(:,kk)));
                end
                inner_dots(jj) = inner_dot;
            end
            
            [~, idxs] = sort(inner_dots);
            new_order(first:last) = idxs(end-degeneracy(num)+1:end);
        end
    end
end

