function vv = prox_quad(vec, inv, AA, bb, eta)
    vv = inv * ( AA'*bb + 1/eta * vec );
    % vv = ( A'*A + 1/eta *eye(length(b)) ) \ ( A'*b + 1/eta * vec ); 
end