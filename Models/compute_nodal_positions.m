function compute_nodal_positions(rod_length, des_face, incl_payload)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%note: 0-indexed below copied from NTRT
Faces = [
    0,8,4;
    0,8,9;
    0,5,9;
    0,3,5;
    0,3,4;
    3,4,11;
    3,10,11;
    3,5,10;
    2,7,11;
    2,10,11;
    2,6,10;
    1,2,6;
    1,2,7;
    1,7,8;
    1,8,9;
    1,6,9;
    4,7,8;
    5,6,9;
    5,6,10;
    4,7,11] + 1;

% construct 6-bar outer structure 
d = rod_length/2;
node_pos = [d, d/2, 0, ...
    -d, d/2, 0, ...
    d, -d/2, 0, ...
    -d, -d/2, 0, ...
    -d/2, 0, -d, ...
    -d/2, 0, d, ...
    d/2, 0, -d, ...
    d/2, 0 , d, ...
    0, -d, d/2, ...
    0, d, d/2, ...
    0, -d, -d/2, ...
    0, d, -d/2
    ]';

% rotate nodes so desired face is bottom
col_node_pos = reshape(node_pos, 3, []);
face_vec = mean([col_node_pos(:,(Faces(des_face,1))), ...
                    col_node_pos(:,(Faces(des_face,2))), ...
                    col_node_pos(:,(Faces(des_face,3)))
                    ], 2);
down_vec = [0, 0, -1]';
% z-y-z
th1 = atan2d(face_vec(2),face_vec(1));
M1z = RRz(-th1);
th2 = atan2d(down_vec(2),down_vec(1));
M2z = RRz(-th2);
v1 = M1z*face_vec;
v2 = M2z*down_vec;
b = atan2d(v2(1),v2(3));
a = atan2d(v1(1),v1(3));
My = RRy(b-a);
R = M2z'*My*M1z;

col_node_pos = R*col_node_pos;
disp('Rotated Node Positions:')
reshape(col_node_pos, [], 1)

end

function R = RRx(theta)
% ccw rotation in DEGREES with rotation axis out of the page
R =  [1  0           0
      0  cosd(theta) -sind(theta)
      0  sind(theta)  cosd(theta)];
end
function R = RRy(theta)
% ccw rotation in DEGREES with rotation axis out of the page
R = [cosd(theta) 0  sind(theta)
     0           1  0              
     -sind(theta) 0  cosd(theta)];
end
function R = RRz(theta)
% ccw rotation in DEGREES with rotation axis out of the page
R =  [cosd(theta) -sind(theta) 0
      sind(theta)  cosd(theta) 0
      0            0           1];
end  


