function [ Rslerped ] = slerpRotationMatrices(R1, R2, w )
%SLERPROTATIONMATRICES Summary
%   Detailed explanation goes here

%% Top-level Algorithm

% convert rotation matrices to axisAngles
[ R1_aa ] = rotation2axisAngle (R1)
[ R2_aa ] = rotation2axisAngle (R2)
% TODO: is there a matlab function for the above to check this?
% vrrotmat2vec(R)?

% convert axisAngles to quaternions
R1_Q = axisAngle2quaternion ( R1_aa )
R2_Q = axisAngle2quaternion ( R2_aa )
% TODO: is there a matlab function for the above to check this?

% slerp between the two quaternions
Q = slerpQuaternions ( R1_Q, R2_Q, w )
% TODO: is there a matlab function for the above to check this?

% convert quaternion to axisAngle
[ aa ] = quaternion2axisAngle ( Q )
% TODO: is there a matlab function for the above to check this?

% convert axisAngle to rotation matrix
Rslerped = axisAngle2rotation ( aa );
% TODO: is there a matlab function for the above to check this?
% vrrotvec2mat([ rn(1), rn(2), rn(3), theta ]) ?
end

function [ aa ] = rotation2axisAngle (R)
%ROTATION2AXISANGLE Summary
%   Detailed explanation goes here

% get eigen vectors and values
[V, D] = eig(R);

% set axis to eigenVector column with corresponding value of 1
[ row, col ] = find(abs(1-D) < 0.0001); % TODO: less wasteful ways to do this? (lamda mod size(D,2))?
rn = V(:, col);

% find theta
cosTheta = (trace(R)-1) / 2;
sincTheta = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]' ./ (2*rn);
sinTheta = sincTheta * norm(rn);

theta = atan2(cosTheta, sinTheta); % TODO: seems to be returning a 3x1 matrix. Should this be a single number?
theta = theta(1); % TODO: brut-forcing for now, numbers seem slightly off as well

aa = [rn; theta];
end

function [ Q ] = axisAngle2quaternion ( aa )
%AXISANGLE2QUATERNION Summary
%   Detailed explanation goes here

theta = aa(4);
sinThetaOver2 = sin(theta / 2);
q1 = aa(1) * sinThetaOver2;
q2 = aa(2) * sinThetaOver2;
q3 = aa(3) * sinThetaOver2;
q4 = cos(aa(4) / 2);

Q = [ q1, q2, q3, q4 ];
end

function [ Qslerped ] = slerpQuaternions ( Q1, Q2, w )
%SLERP Summary
%   Detailed explanation goes here

phi = acos(Q1 .* Q2);
Qslerped = (Q1 * sin(1-w) .* phi./sin(phi)) + (Q2.*sin(w * phi)./sin(phi));
end

function [ aa ] = quaternion2axisAngle ( Q )
%QUATERNION2AXISANGLE Summary
%   Detailed explanation goes here
theta = 2 * acos(Q(4))
sinThetaOver2 = sin(theta / 2)
rn = [                          ...
    Q(1) ./ sinThetaOver2    ;   ...
    Q(2) ./ sinThetaOver2    ;   ...
    Q(3) ./ sinThetaOver2        ...
    ]

aa = [ rn; theta ]
end

function [ R ] = axisAngle2rotation ( aa )
%AXISANGLE2ROTATION Summary
%   Detailed explanation goes here
theta = aa(4);
rn = [ aa(1); aa(2); aa(3) ];
r = rn * theta;
rMag = norm(r); % TODO: doesn't this just equal theta?
rx = [                          ...
    0       -r(3)	r(2)    ;   ...
    r(3)    0       -r(1)   ;   ...
    -r(2)   r(1)    0           ...
    ];

R = cos(rMag) * eye(3)                  ...
    + sin(rMag) / rMag .* rx            ...
    + ((1-cos(rMag)) / rMag^2) * r * r' ;
end
