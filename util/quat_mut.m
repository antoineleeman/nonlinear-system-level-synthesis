function ret = quat_mut(p,q)
    p0 = p(1);
    p = p(2:4);
    q0 = q(1);
    q = q(2:4);
    ret = [p0*q0 - p'*q;
        p0*q + q0*p + cross(p,q)];
end