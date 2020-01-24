%Create the graph from the point cloud using an epsilon-neighborhood
%connectivity
clear
bunnyclean = gsp_pointcloud('bunny');
noise = wgn(2503,1,-60);
noise = repmat(noise,1,3);
bunnynoise = bunnyclean + noise;
param.type = 'knn';
param.rescale = 1;
param.center = 1;

%Compute it
G1 = gsp_nn_graph(double(bunnyclean), param);
%Reduce vertex size for plotting
G1.plotting.vertex_size = 10;

%Compute it
G2 = gsp_nn_graph(double(bunnynoise), param);
%Reduce vertex size for plotting
G2.plotting.vertex_size = 10;
% G1 = gsp_bunny;
% G2 = gsp_bunny;
figure;
gsp_plot_graph(G1)
figure;
gsp_plot_graph(G2)
%G1 = gsp_random_regular();
%G2 = gsp_random_regular();