function vv = prox_l1(ww, ll)
    I = abs(ww) > ll;
    vv = zeros(size(ww));

    vv(I) = ww(I) - ll* sign(ww(I));
end