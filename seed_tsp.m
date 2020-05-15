

%
% Create XY points
%
xy = 10*rand(200,2);


%
% Solve the TSP using NN
%
nn = tsp_nn('xy',xy);


%
% Solve the TSP starting with the NN solution
%
ga = tsp_ga(nn);


%
% Compare against running the TSP from scratch
%
ga0 = tsp_ga('xy',xy);

