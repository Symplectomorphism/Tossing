function w = skewvec(R)

w = zeros(3,1);
w(1) = -R(2,3);
w(2) = R(1,3);
w(3) = -R(1,2);