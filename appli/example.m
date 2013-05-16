%EXAMPLE

% input is a 72x12 feature matrix.
load 'input.mat';

% First 36 rows are of class A, last 36 rows are of class B
target=[ones(36,1) ; 2*ones(36,1)];

% change mode from 1 to 0 to see what happens.
matlab2csv(input,target,2);

