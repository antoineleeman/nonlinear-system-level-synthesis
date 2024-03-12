function [q] = euler_to_quat(a)
a = deg2rad(a);
% pitch, roll, yaw
cy = cos(a(2) * 0.5);
sy = sin(a(2) * 0.5);
cr = cos(a(1) * 0.5);
sr = sin(a(1) * 0.5);
cp = cos(a(3) * 0.5);
sp = sin(a(3) * 0.5);

q= zeros(4,1);

q(1) = cy * cr * cp + sy * sr * sp;
q(2) = cy * sr * cp - sy * cr * sp;
q(3) = cy * cr * sp + sy * sr * cp;
q(4) = sy * cr * cp - cy * sr * sp;

end

