function hp = pcolorx(e1, e2, vals)

de1 = e1(2) - e1(1);
de2 = e2(2) - e2(1);

xe1 = reshape(e1, numel(e1), 1);
xe2 = reshape(e2, numel(e2), 1);

xe1 = [xe1 - de1/2; xe1(end) + de1/2];
xe2 = [xe2 - de2/2; xe2(end) + de2/2];

xvals = vals;
xvals = [xvals; xvals(1,:)];
xvals = [xvals, xvals(:,1)];

hp = pcolor(xe1, xe2, xvals);
