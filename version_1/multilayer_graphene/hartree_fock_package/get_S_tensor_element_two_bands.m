function S_element = get_S_tensor_element_two_bands(eigvec1, valley1, eigvec2, valley2, eigvec3, valley3, eigvec4, valley4, LL_index_cutoff, d_interlayer)
    % evaluate the S matrix @ (\lambda, \sigma ; \alpha, \beta) for bilayer graphene
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，维数为dims = 2 * LL_index_cutoff
    % valley : + 1 for K valley ; -1 for Kp valley
    % in K valley, the Landau levels in the order : [|0,B2>, |1,B2> // |2,B2>, |0,A1> // |3,B2>, |1,A1> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,A1>, |1,A1> // |2,A1>, |0,B2> // |3,A1>, |1,B2> // ... ...]
    
    S_element = 0.0;
    if ~(valley1 == valley4) % 在这里我们不考虑太复杂的intervalley term, 因为它们发生的概率很小
        return
    end
    
    if ~(valley2 == valley3) % 在这里我们不考虑太复杂的intervalley term, 因为它们发生的概率很小
        return
    end
    
    [m1_more_list, m1_less_list] = get_LL_index_each_sublattice(LL_index_cutoff, valley1);
    m4_more_list = m1_more_list;
    m4_less_list = m1_less_list;
    m1_index_bound = length(m1_less_list);
    
    [m2_more_list, m2_less_list] = get_LL_index_each_sublattice(LL_index_cutoff, valley2);
    m3_more_list = m2_more_list;
    m3_less_list = m2_less_list;
    m2_index_bound = length(m2_less_list);
    
    % less & less
    for index1 = 1:m1_index_bound
        for index4 = 1:m1_index_bound
            for index2 = 1:m2_index_bound
                for index3 = 1:m2_index_bound
                    m1 = index1 - 1;
                    m4 = index1 - 1;
                    m2 = index1 - 1;
                    m3 = index1 - 1;
                    
                    coeff_intralayer = conj(eigvec1(m1_more_list(index1))) * eigvec4(m4_more_list(index4)) * conj(eigvec2(m2_more_list(index2))) * eigvec3(m3_more_list(index3));
                    coeff_intralayer = coeff_intralayer + conj(eigvec1(m1_less_list(index1))) * eigvec4(m4_less_list(index4)) * conj(eigvec2(m2_less_list(index2))) * eigvec3(m3_less_list(index3));
                    coeff_interlayer = conj(eigvec1(m1_more_list(index1))) * eigvec4(m4_more_list(index4)) * conj(eigvec2(m2_less_list(index2))) * eigvec3(m3_less_list(index3));
                    coeff_interlayer = coeff_interlayer + conj(eigvec1(m1_less_list(index1))) * eigvec4(m4_less_list(index4)) * conj(eigvec2(m2_more_list(index2))) * eigvec3(m3_more_list(index3));
                    
                    X_nm_intralayer = get_exchange_integral(m1, m4, m2, m3); % intralayer
                    X_nm_interlayer = get_exchange_integral(m1, m4, m2, m3, d_interlayer); % interlayer
                    
                    S_element = S_element + coeff_intralayer * X_nm_intralayer;
                    S_element = S_element + coeff_interlayer * X_nm_interlayer;
                end
            end
        end
    end
    
    % more & more
    if valley1 == valley2
        for index1 = m1_index_bound + 1:m1_index_bound + 2
            for index4 = m1_index_bound + 1:m1_index_bound + 2
                for index2 = m2_index_bound + 1:m2_index_bound + 2
                    for index3 = m2_index_bound + 1:m2_index_bound + 2
                        m1 = index1 - 1;
                        m4 = index1 - 1;
                        m2 = index1 - 1;
                        m3 = index1 - 1;
                        
                        coeff_intralayer = conj(eigvec1(m1_more_list(index1))) * eigvec4(m4_more_list(index4)) * conj(eigvec2(m2_more_list(index2))) * eigvec3(m3_more_list(index3));
                        coeff_intralayer = coeff_intralayer + conj(eigvec1(m1_less_list(index1))) * eigvec4(m4_less_list(index4)) * conj(eigvec2(m2_less_list(index2))) * eigvec3(m3_less_list(index3));
                        
                        X_nm_intralayer = get_exchange_integral(m1, m4, m2, m3); % intralayer
                        
                        S_element = S_element + coeff_intralayer * X_nm_intralayer;
                    end
                end
            end
        end
    else
        for index1 = m1_index_bound + 1:m1_index_bound + 2
            for index4 = m1_index_bound + 1:m1_index_bound + 2
                for index2 = m2_index_bound + 1:m2_index_bound + 2
                    for index3 = m2_index_bound + 1:m2_index_bound + 2
                        m1 = index1 - 1;
                        m4 = index1 - 1;
                        m2 = index1 - 1;
                        m3 = index1 - 1;
                        
                        coeff_interlayer = conj(eigvec1(m1_more_list(index1))) * eigvec4(m4_more_list(index4)) * conj(eigvec2(m2_less_list(index2))) * eigvec3(m3_less_list(index3));
                        coeff_interlayer = coeff_interlayer + conj(eigvec1(m1_less_list(index1))) * eigvec4(m4_less_list(index4)) * conj(eigvec2(m2_more_list(index2))) * eigvec3(m3_more_list(index3));
                    
                        X_nm_interlayer = get_exchange_integral(m1, m4, m2, m3, d_interlayer); % interlayer
                        
                        S_element = S_element + coeff_interlayer * X_nm_interlayer;
                    end
                end
            end
        end
    end
    
    % more & less
    if valley1 == valley2
        
    else
        
    end
end


function [m_more_list, m_less_list] = get_LL_index_each_sublattice(LL_index_cutoff, valley)
    if valley == 1 % for K valley
        m_less_list = zeros(1, LL_index_cutoff - 1);
        m_more_list = zeros(1, 2 + LL_index_cutoff - 1);
        m_more_list(1) = 1;
        m_more_list(2) = 2;
        
        for n = 2:LL_index_cutoff
            m_more_list(n + 1) = 2 * (n - 1) + 1;
            m_less_list(n - 1) = 2 * n;
        end
   
    else % for Kp valley
        m_more_list = zeros(1, LL_index_cutoff + 1);
        m_less_list = zeros(1, LL_index_cutoff - 1);
        m_more_list(1) = 1;
        m_more_list(2) = 2;
        
        for n = 2:LL_index_cutoff
            m_more_list(n + 1) = 2 * (n - 1) + 1;
            m_less_list(n - 1) = 2 * n;
        end
    end
end