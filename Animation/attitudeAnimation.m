function [] = attitudeAnimation(A1, A2, stlName, percent_start, step_size, animation_length)

% creates an animation with 2 attitudes
%
% [] = attitudeAnimation(A1, A2, stlName, percent_start, step_size, animation_lenght)
%
% Input arguments:
% ----------------------------------------------------------------
% A1                    [3x3xN]         A matrix n1             [-]
% A2                    [3x3xN]         A matrix n2             [-]
% stlName               [1x1]           Location of .stl file   [str]
% percent_start         [1x1]           Where to start          [%]
% step_size             [1x1]           Step size               [-]
% animation_length      [1x1]           Total states computed   [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% N/A

[TR] = stlImport_internal(stlName);

if length(A1) ~= length(A2)
    error("A1 and A2 must have the same length!")
end

A_len = length(A1);
if percent_start == 0
    n_start = 1;
else
    n_start = ceil(A_len * percent_start * 1e-2);
end
if animation_length == 0
    n_end = A_len;
else
    n_end = n_start + animation_length;
    if n_end > A_len
        n_end = A_len;
    end
end

% calculate rotations
k = 0;
TR_frame_vect1(1:length(n_start:step_size:n_end), 1) = {NaN(1,1)};
TR_frame_vect2(1:length(n_start:step_size:n_end), 1) = {NaN(1,1)};
for i = n_start:step_size:n_end
    k = k + 1;
    [TR_frame1] = stlRotate_internal(TR, A1(1:3,1:3,i));
    TR_frame_vect1{k} = TR_frame1;
    [TR_frame2] = stlRotate_internal(TR, A2(1:3,1:3,i));
    TR_frame_vect2{k} = TR_frame2;
end

% calculate max frame
x_min = min(TR_frame_vect1{1}.Points(:,1));
x_max = max(TR_frame_vect1{1}.Points(:,1));
y_min = min(TR_frame_vect1{1}.Points(:,2));
y_max = max(TR_frame_vect1{1}.Points(:,2));
z_min = min(TR_frame_vect1{1}.Points(:,3));
z_max = max(TR_frame_vect1{1}.Points(:,3));
for i = 1:length(n_start:step_size:n_end)
    x_min_1 = min(TR_frame_vect1{i}.Points(:,1));
    x_max_1 = max(TR_frame_vect1{i}.Points(:,1));
    y_min_1 = min(TR_frame_vect1{i}.Points(:,2));
    y_max_1 = max(TR_frame_vect1{i}.Points(:,2));
    z_min_1 = min(TR_frame_vect1{i}.Points(:,3));
    z_max_1 = max(TR_frame_vect1{i}.Points(:,3));
    x_min_2 = min(TR_frame_vect2{i}.Points(:,1));
    x_max_2 = max(TR_frame_vect2{i}.Points(:,1));
    y_min_2 = min(TR_frame_vect2{i}.Points(:,2));
    y_max_2 = max(TR_frame_vect2{i}.Points(:,2));
    z_min_2 = min(TR_frame_vect2{i}.Points(:,3));
    z_max_2 = max(TR_frame_vect2{i}.Points(:,3));
    if min(x_min_1, x_min_2) < x_min
        x_min = min(x_min_1, x_min_2);
    end
    if max(x_max_1, x_max_2) > x_max
        x_max = max(x_max_1, x_max_2);
    end
    if min(y_min_1, y_min_2) < y_min
        y_min = min(y_min_1, y_min_2);
    end
    if max(y_max_1, y_max_2) > y_max
        y_max = max(y_max_1, y_max_2);
    end
    if min(z_min_1, z_min_2) < z_min
        z_min = min(z_min_1, z_min_2);
    end
    if max(z_max_1, z_max_2) > z_max
        z_max = max(z_max_1, z_max_2);
    end
end

h = figure("WindowState", "maximized");

start = 1;
k = 0;
for i=n_start:step_size:n_end
    k = k + 1;

    subplot(1,2,1)
    title("Attitude 1")
    xlim([x_min x_max])
    ylim([y_min y_max])
    zlim([z_min z_max])

    a1 = trisurf(TR_frame_vect1{k}, 'FaceColor','black','EdgeColor','interp');

    if start
        hold on, axis equal
        % light
    end

    vect = 1500*[A1(1:3,1,i)'; 0,0,0];
    b1 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="red"); % x
    vect = 1500*[A1(1:3,2,i)'; 0,0,0];
    c1 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="green"); % y
    vect = 1500*[A1(1:3,3,i)'; 0,0,0];
    d1 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="blu"); % z

    subplot(1,2,2)
    title("Attitude 2")
    xlim([x_min x_max])
    ylim([y_min y_max])
    zlim([z_min z_max])

    a2 = trisurf(TR_frame_vect2{k}, 'FaceColor','black','EdgeColor','interp');
    
    if start
        hold on, axis equal
        start = 0;
        % light
    end
    
    vect = 1500*[A2(1:3,1,i)'; 0,0,0];
    b2 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="red"); % x
    vect = 1500*[A2(1:3,2,i)'; 0,0,0];
    c2 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="green"); % y
    vect = 1500*[A2(1:3,3,i)'; 0,0,0];
    d2 = plot3(vect(:,1), vect(:,2), vect(:,3), Color="blu"); % z
    
    drawnow
    pause(0.05)
    delete(a1)
    delete(b1)
    delete(c1)
    delete(d1)
    delete(a2)
    delete(b2)
    delete(c2)
    delete(d2)
end

close(h)

end


%% internal functions

function [TR] = stlImport_internal(stl_file_name)
    gm = stlread(stl_file_name);
    
    x_dim = max(gm.Points(:,1))-min(gm.Points(:,1));
    y_dim = max(gm.Points(:,2))-min(gm.Points(:,2));
    z_dim = max(gm.Points(:,3))-min(gm.Points(:,3));
    
    % centre of the body
    cm_body = [x_dim/2 y_dim/2 z_dim/2];
    transposition = [min(gm.Points(:,1)) min(gm.Points(:,2)) min(gm.Points(:,3))];
    
    % points with centre of the body in [0 0 0]
    new_points = gm.Points - ones(size(gm.Points)).*(cm_body + transposition);
    TR = triangulation(gm.ConnectivityList, new_points);
end
function [TR] = stlRotate_internal(TR, Rot_matrix)
    newTrPoints = zeros(size(TR.Points));
    for i = 1:size(TR.Points, 1)
        newTrPoints(i,:) = Rot_matrix * TR.Points(i,:)';
    end
    TR = triangulation(TR.ConnectivityList, newTrPoints);
end

