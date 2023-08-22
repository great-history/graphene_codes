function H_hf_new = refine_H_hf(H_hf, dim)
    % 虽然理论上是厄密的，但在构造full hartree fock ham的时候，可能会有一些偏差(比如对角元上出现一个很小的虚部)
    % 因此重新使其严格厄米
    H_hf_new = zeros(size(H_hf));
    for ii = 1:dim
        for jj = ii:dim
            H_hf_new(ii, jj) = (H_hf(ii, jj) + conj(H_hf(jj, ii))) / 2;
            H_hf_new(jj, ii) = conj(H_hf_new(ii, jj));
        end
    end
    
end