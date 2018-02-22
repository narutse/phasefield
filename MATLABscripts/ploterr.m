%% Spatial Error:

M = csvread("data/FBEspaceL2.csv");
dx2 = M(:,1);
err = M(:, 2);

X = [ones(length(dx2), 1) log(dx2)];
b = X\log(err);
title('log log plot of spatial error')
loglog(dx2, err, '-o')

%% Time Error

M = csvread("data/FBEtimeL2.csv");
dt = M(:, 1);
err =  M(:,2);

X = [ones(length(dt),1) log(dt)];
b = X\log(err);

title('log log plot of temporal error')
loglog(dt, err, '-o')

Xfixed = X(1:8, :);
bfixed = Xfixed\log(err(1:8,:));