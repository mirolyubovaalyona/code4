%=======================Падение по Z
n = 101
m = 0:1:n-1;
theta_ = (pi/(n-1)*m);
X_2 = 0:1:double(n)-1;
fi_ = (pi/2/double(n-1)*m);

theta = repmat(theta_', 1, n);
fi = repmat(fi_, n, 1);

F = 1

X = F .* sin(theta) .* cos(fi);
Y = F .* sin(theta) .* sin(fi);
Z = F .* cos(theta);

figure
hold on;

surf(X, Y, Z);
surf(X, -Y, Z);
surf(-X, Y, Z);
surf(-X, -Y, Z);
axis square
title('Падение по Z')
hold off;

%=======================Падение по X

n = 101
m = 0:1:n-1;
theta_ = (pi/2/(n-1)*m);
fi_ = (pi/double(n-1)*m);

theta = repmat(theta_', 1, n);
fi = repmat(fi_, n, 1);

F = 1

X= F .* sin(theta) .* cos(fi);
Y = F .* sin(theta) .* sin(fi);
Z = F .* cos(theta);

figure
hold on;

surf(X, Y, Z);
surf(X, -Y, Z);
surf(X, Y, -Z);
surf(X, -Y, -Z);
axis square
title('Падение по X')
hold off;

%======================= Падение по Y 

n = 101
m = 0:1:n-1;
theta_ = (pi/2/(n-1)*m);
X_2 = 0:1:double(n)-1;
fi_ = (pi/2 + pi/double(n-1)*m);

theta = repmat(theta_', 1, n);
fi = repmat(fi_, n, 1);

F = 1

X= F .* sin(theta) .* cos(fi);
Y = F .* sin(theta) .* sin(fi);
Z = F .* cos(theta);

figure
hold on;
axis square
surf(X, Y, Z);
surf(-X, Y, Z);
surf(X, Y, -Z);
surf(-X, Y, -Z);

title('Падение по Y')
hold off;