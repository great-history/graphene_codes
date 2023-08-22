% 判断一个矩阵的矩阵元是否都是实数

function hem = helper_check_real(matrix, eps)
    m_size = size(matrix);
    hem = 1;
    for i = 1:m_size(1)
        for j = 1:m_size(2)
            if abs(imag(matrix(i,j))) >= eps
                hem = 0;
            end
        end
    end
end