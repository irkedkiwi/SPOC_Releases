% pdf plot example
% For a 2d problem...


Discretization = [50, 200,200];
Methods = [1,0,0];
[time_nodes, parameter_values, pdf_values] = ... 
                                    Extract_Nodes(Discretization, Methods);
    
cmap=flipud(hot(256));


f_colors = cell(length(parameter_values(:,1)),2);
for mesh_index = 1:length(length(parameter_values(:,1)))
    yy = linspace(0,1,size(cmap,1));    % Generate range of color indices that map to cmap
    cm = spline(yy,cmap',pdf_values(:));                  % Find interpolated colorvalues 
    cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
    cm(cm<0)=0;
end
 
figure
hold on
for mesh_index = 1:length(parameter_values(:,1))
    plot(parameter_values(mesh_index,1),parameter_values(mesh_index,2),'.','color',cm(:,mesh_index),'LineWidth',1.5)
end

