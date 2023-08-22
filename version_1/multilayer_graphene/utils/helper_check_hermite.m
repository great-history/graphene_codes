% 判断一个矩阵是否是厄密矩阵

function hem = helper_check_hermite(matrix, eps)
    m_size = size(matrix);
    hem = 1;
    if m_size(1) == m_size(2)
        for i = 1:(m_size(1)-1)
            for j = (i+1):m_size(1)
                error = abs(matrix(i,j) - conj(matrix(j,i)));
                if error > eps
                    hem = 0;
                end
            end
        end
    else
        hem = 0;
    end
end