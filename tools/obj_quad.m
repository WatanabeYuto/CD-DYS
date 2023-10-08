function f = obj_quad(A,b,x)
    f = 0.5 * (A*x - b)'*(A*x - b);
end