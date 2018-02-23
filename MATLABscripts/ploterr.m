%% Spatial Error:

M = csvread("../data/FBEspaceL2.csv");
dx2 = M(:,1);
err = M(:, 2);

X = [ones(length(dx2), 1) log(dx2)];
b_sp = X\log(err);
loglog(dx2, err, '-o')
title('Spatial Error (loglog)')

%% Time Error

M = csvread("../data/FBEtimeL2.csv");
dt = M(:, 1);
err =  M(:,2);


loglog(dt, err, '-o')
title('Temporal error (loglog)')

X = [ones(length(dt),1) log(dt)];
b_t = X\log(err);

Xfixed = X(1:8, :);
bfixed = Xfixed\log(err(1:8,:));