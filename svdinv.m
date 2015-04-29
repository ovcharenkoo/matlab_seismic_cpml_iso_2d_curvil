%SVD-like matrix inversion
function B=svdinv(A)
    [U, S, V]= svd(A);
    s= diag(S); 
    k= sum(s> 1e-9); % simple thresholding based decision
    B= V(:, 1: k)* diag(1./ s(1: k))* U(:, 1: k)';
end