%% the input of the beam element solver
path = 'D:\OneDrive - McGill University\publication\Journal of mechanical science\bending\cubic_input\';
fileID = fopen(strcat(path,'topnode.txt'),'r');
formatSpec = '%d';
top = fscanf(fileID,formatSpec);
top = top';
fileID = fopen(strcat(path,'bottomnode.txt'),'r');
formatSpec = '%d';
bottom = fscanf(fileID,formatSpec);
bottom = bottom';
fixdofs = [bottom*6-5 bottom*6-4 bottom*6-3 bottom*6-2 bottom*6-1 bottom*6 ...
    top*6-5 top*6-4 top*6-3 top*6-2 top*6-1 top*6];
load = zeros(length(top),2);
load(:,1) = top*6-3;
load(:,2) = -1000;
fileID = fopen(strcat(path,'node_list.txt'),'r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
node_list = fscanf(fileID,formatSpec,sizeA);
node_list = node_list';
fileID = fopen(strcat(path,'strut_list.txt'),'r');
formatSpec = '%d %d';
sizeA = [2 Inf];
strut_list = fscanf(fileID,formatSpec,sizeA);
strut_list = strut_list';






