function ginv = homInv(g)

R = g(1:3,1:3);
p = g(1:3, 4);

ginv = [R', -(R')*p; zeros(1,3), 1];

end