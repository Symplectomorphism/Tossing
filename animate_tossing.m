function animate_tossing(t, q)

% Screen bounds
xmin = min( q(:,1) ) - 2 * max( q(:,25) );
ymin = min( q(:,2) ) - 2 * max( q(:,25) );
zmin = min( q(:,3) ) - 2 * max( q(:,25) );
xmax = max( q(:,1) ) + 2 * max( q(:,25) );
ymax = max( q(:,2) ) + 2 * max( q(:,25) );
zmax = max( q(:,3) ) + 2 * max( q(:,25) );

lent = length(t);
stepSize = 1;

scrsz = get(groot, 'ScreenSize');
fig = figure(100);
clf(fig);
set(fig, 'OuterPosition', [1 scrsz(4)/6 5*scrsz(3)/6 5*scrsz(4)/6]);

etaspan = 0:pi/72:2*pi;
dough = zeros(length(etaspan),3);
i = 1;
for eta = etaspan
    dough(i,:) = [cos(eta), sin(eta), 0];
    i = i + 1;
end


hd = fill3(dough(:,1),dough(:,2),dough(:,3), [0.6,0.6,0.6]);
hl = line([0, q(1,25)],[0,0],[0,0], 'Color', 'red', 'LineWidth', 2);
% hs = patch(1,1,1, 'FaceColor', 0.6*ones(3,1));

axis([xmin, xmax, ymin, ymax, zmin, zmax])

%%%-----------------------------------------------------------------%%%
% Draw Workspace
%----------------------------------------------------------------------
% Bottom
line( [xmin; xmax], [ymin; ymin], [zmin, zmin], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Top
line( [xmin; xmax], [ymax; ymax], [zmin, zmin], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Left
line( [xmin; xmin], [ymin; ymax], [zmin, zmin], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Right
line( [xmax; xmax], [ymin; ymax], [zmin, zmin], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')

% Bottom
% line( [xmin; xmin], [ymin; ymin], [zmin, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Top
line( [xmin; xmin], [ymax; ymax], [zmin, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Left
line( [xmax; xmax], [ymax; ymax], [zmin, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% Right
line( [xmin; xmax], [ymax; ymax], [zmax, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
% line( [xmin; xmin], [ymin; ymax], [zmax, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')

line( [xmax; xmax], [ymin; ymin], [zmin, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
line( [xmax; xmax], [ymin; ymax], [zmax, zmax], 'LineWidth', 2, 'Color', [.1 .1 .1], 'LineStyle', '--')
%----------------------------------------------------------------------

xlabel('x')
ylabel('y')
zlabel('z')

for i = 1:stepSize:lent
    try
        delete(hs)
    end
    R = reshape(q(i, 7:15), 3, 3)';
    
    
    temp = q(i,1:3)' + R*[q(i,25); 0; 0];
    hl.XData = [q(i,1), temp(1)];
    hl.YData = [q(i,2), temp(2)];
    hl.ZData = [q(i,3), temp(3)];
    
    if abs(q(i,27)) > 1e-6
        delete(hd)
        
        xi = linspace(0,q(i,27),21);
        eta = 0:pi/18:2*pi;
        [xi, eta] = meshgrid(xi, eta);
        xt = q(i,25)*sqrt(abs(q(i,27)-xi)/abs(q(i,27))).*cos(eta);
        yt = q(i,25)*sqrt(abs(q(i,27)-xi)/abs(q(i,27))).*sin(eta);
        zt = xi;
        
        xd = zeros(size(xt));
        yd = zeros(size(yt));
        zd = zeros(size(zt));
        for j = 1:size(xt,1)
            for k = 1:size(xt,2)
                temp = q(i,1:3)' + R*[xt(j,k); yt(j,k); zt(j,k)];
                xd(j,k) = temp(1);
                yd(j,k) = temp(2);
                zd(j,k) = temp(3);
            end
        end
        fvc = surf2patch(xd,yd,zd);
        hs = patch(fvc, 'FaceColor', [0.3, 0.3, 0.3], 'EdgeColor', 'g');
    else
        delete(hd);
        hd = patch(dough(:,1),dough(:,2),dough(:,3), [0.6,0.6,0.6]);
        j = 1;
        temp = zeros(length(etaspan),3);
        for eta = etaspan
            temp(j,:) = (R*q(i,25)*dough(j,:)')';
            j = j + 1;
        end
        hd.XData = q(i,1) + temp(:,1)';
        hd.YData = q(i,2) + temp(:,2)';
        hd.ZData = q(i,3) + temp(:,3)';
    end
    
%     pause(0.05)
    drawnow
    view(-45, 45)
    lighting gouraud
end

end