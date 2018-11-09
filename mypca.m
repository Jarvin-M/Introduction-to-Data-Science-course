function [pc, eigenvalues] = mypca(A)
%A is a m X n matrix with m observations and n variables.
%pc represents the L2-normalized principal components
%eigenvalues represent the corresponding eigen values of the principal
%components

    %Center the matrix by subtracting the mean from each component
    centeredA = A-mean(A);
    %Find the covariance of the centered values of A
    covarianceMatrix = cov(centeredA);
    
    %function eig returns the eigenvalues and eigenvectors of the covariance of the
    %centered matrix
    [eigenvectors,~] = eig(covarianceMatrix);
    eigenvalues = eig(covarianceMatrix);
    
    %pc - L2-Normalized components 
    pc = zeros(size(eigenvectors));
    for i = 1:size(eigenvectors,2)
        
        pc(:,i)= eigenvectors(:,i)/sqrt(sum(eigenvectors(:,i).*eigenvectors(:,i)));
    end
    
    %sorting
    [~, srtidx] = sort(eigenvalues, 'descend');
    eigenvalues = sort(eigenvalues, 'descend');
    pc = pc(:,srtidx); 
end

