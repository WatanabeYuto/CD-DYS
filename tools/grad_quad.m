function g = grad_quad(A,b,x)
    %% the gradient of \frac{1}{2} \|Ax - b\|^2
    g = A'*(A*x-b);
end