function [] = attitudeAnimation2(A, P, theta, stlName, percent_start, step_size, animation_length, cam_choice)

% creates an animation with attitude and position
%
% [] = attitudeAnimation(A1, A2, stlName, percent_start, step_size, animation_lenght)
%
% Input arguments:
% ----------------------------------------------------------------
% A1                    [3x3xN]     A matrix                        [-]
% P                     [Nx3]       Position matrix                 [-]
% theta                 [Nx1]       Theta angle matrix              [-]
% stlName               [1x1]       Location of .stl file           [str]
% percent_start         [1x1]       Where to start                  [%]
% step_size             [1x1]       Step size                       [-]
% animation_length      [1x1]       Total states computed           [-]
% cam_choice            [1x1]       1 = facing Earth, 2 = facing h  [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% N/A

warning("Earth is not rotating")

[TR] = stlImport_internal(stlName);

A_len = length(A);
if percent_start == 0
    n_start = 1;
else
    n_start = ceil(A_len * percent_start * 1e-2);
end
n_end = n_start + animation_length;
if n_end > A_len
    n_end = A_len;
end

% calculate rotations
k = 0;
TR_frame_vect(1:length(n_start:step_size:n_end), 1) = {NaN(1,1)};
for i = n_start:step_size:n_end
    k = k + 1;
    [TR_frame1] = stlHandle_internal(TR, A(1:3,1:3,i), P(i,:));
    TR_frame_vect{k} = TR_frame1;
end

% calculate max frame
x_min = min(TR_frame_vect{1}.Points(:,1));
x_max = max(TR_frame_vect{1}.Points(:,1));
y_min = min(TR_frame_vect{1}.Points(:,2));
y_max = max(TR_frame_vect{1}.Points(:,2));
z_min = min(TR_frame_vect{1}.Points(:,3));
z_max = max(TR_frame_vect{1}.Points(:,3));
for i = 1:length(n_start:step_size:n_end)
    x_min_1 = min(TR_frame_vect{i}.Points(:,1));
    x_max_1 = max(TR_frame_vect{i}.Points(:,1));
    y_min_1 = min(TR_frame_vect{i}.Points(:,2));
    y_max_1 = max(TR_frame_vect{i}.Points(:,2));
    z_min_1 = min(TR_frame_vect{i}.Points(:,3));
    z_max_1 = max(TR_frame_vect{i}.Points(:,3));
    if x_min_1 < x_min
        x_min = x_min_1;
    end
    if x_max_1 > x_max
        x_max = x_max_1;
    end
    if y_min_1 < y_min
        y_min = y_min_1;
    end
    if y_max_1 > y_max
        y_max = y_max_1;
    end
    if z_min_1 < z_min
        z_min = z_min_1;
    end
    if z_max_1 > z_max
        z_max = z_max_1;
    end
end
earth_x_min = -6371.01;
earth_x_max = 6371.01;
earth_y_min = -6371.01;
earth_y_max = 6371.01;
earth_z_min = -6371.01;
earth_z_max = 6371.01;
if x_min > earth_x_min
    x_min = earth_x_min;
end
if x_max < earth_x_max
    x_max = earth_x_max;
end
if y_min > earth_y_min
    y_min = earth_y_min;
end
if y_max < earth_y_max
    y_max = earth_y_max;
end
if z_min > earth_z_min
    z_min = earth_z_min;
end
if z_max < earth_z_max
    z_max = earth_z_max;
end

% animation
h = figure("WindowState", "maximized");
title("Attitude")
earthPlot_internal;

h_orbit = cross(P(1,:), P(2,:));
h_orbit = h_orbit./norm(h_orbit);

start = 1;
k = 0;
for i=n_start:step_size:n_end
    k = k + 1;
    a = trisurf(TR_frame_vect{k}, 'FaceColor','black','EdgeColor','interp');

    if start
        hold on, axis equal, axis off
        camzoom(15)
        start = 0;
        % light
    end

    vect = 200*[A(1:3,1,i)'; 0,0,0] + [P(i,:); P(i,:)];
    b = plot3(vect(:,1), vect(:,2), vect(:,3), Color="red"); % x
    vect = 200*[A(1:3,2,i)'; 0,0,0] + [P(i,:); P(i,:)];
    c = plot3(vect(:,1), vect(:,2), vect(:,3), Color="green"); % y
    vect = 200*[A(1:3,3,i)'; 0,0,0] + [P(i,:); P(i,:)];
    d = plot3(vect(:,1), vect(:,2), vect(:,3), Color="blu"); % z

    xlim([x_min x_max])
    ylim([y_min y_max])
    zlim([z_min z_max])

    if cam_choice == 1
        view(P(i,:))
    elseif cam_choice == 2
        view(h_orbit)
    end
    angle=360*((theta(i+1)-theta(i))/(2*pi));
    camtarget(P(i,:));
    camorbit(angle, 0, 'data', h_orbit)

    drawnow

    pause(0.05)

    delete(a)
    delete(b)
    delete(c)
    delete(d)
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
    new_points = (gm.Points - ones(size(gm.Points)).*(cm_body + transposition)) * 5e-2;     % should be * 1e-3 (was in mm)
    TR = triangulation(gm.ConnectivityList, new_points);
end
function [TR] = stlHandle_internal(TR, Rot_matrix, newPos)
    newTrPoints = zeros(size(TR.Points));
    for i = 1:size(TR.Points, 1)
        newTrPoints(i,:) = Rot_matrix * TR.Points(i,:)';
    end
    TR = triangulation(TR.ConnectivityList, newTrPoints);
    
    % transposition
    x_dim = max(TR.Points(:,1))-min(TR.Points(:,1));
    y_dim = max(TR.Points(:,2))-min(TR.Points(:,2));
    z_dim = max(TR.Points(:,3))-min(TR.Points(:,3));
    
    % reset to position 0
    cm_body = [x_dim/2 y_dim/2 z_dim/2];
    transposition = [min(TR.Points(:,1)) min(TR.Points(:,2)) min(TR.Points(:,3))];
    new_points = TR.Points - ones(size(TR.Points)).*(cm_body + transposition);
    
    new_points = new_points + ones(size(TR.Points)).*(newPos);
    TR = triangulation(TR.ConnectivityList, new_points);
end
function s = earthPlot_internal
    % Earth settings
    [x1,y1,z1] = sphere(200);
    mult=6378;                     % Earth radius
    s = surface(x1*mult,y1*mult,z1*mult);
    
    load topo
    s.FaceColor = 'texturemap';
    s.CData = topo;
    s.EdgeColor = 'none';          % remove edges
    s.FaceLighting = 'gouraud';    % lighting for specific surfaces
    s.SpecularStrength = 0.4;      % reflected light
    axis equal
    
    image_file = "flat_earth.jpg";
    cdata = imread(image_file);
    cdata=cdata(end:-1:1,:,:);
    set(s, 'FaceColor', 'texturemap', 'CData', cdata, 'EdgeColor', 'none');
    
    hold on
    axis equal, grid on
end

