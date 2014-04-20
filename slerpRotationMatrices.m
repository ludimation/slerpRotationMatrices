function [ Rslerped_dallen ] = slerpRotationMatrices( R1, R2, w )
%SLERPROTATIONMATRICES Summary
%   Detailed explanation goes here
%      TODO: might find some more ideas here ? http://www.mathworks.com/matlabcentral/newsreader/view_thread/301137

%% convert rotation matrices to axisAngles
R1_aa_dallen = rotmat2axisAngle ( R1 )
% R1_aa_vrrotmat2vec = vrrotmat2vec( R1 ) % NOTE: matches this MATLAB function value
R2_aa_dallen = rotmat2axisAngle ( R2 )
% R2_aa_vrrotmat2vec = vrrotmat2vec( R2 ) % NOTE: matches this MATLAB function value

%% convert axisAngles to quaternions
%      NOTES: calculations seem correct, but signs are reversed which is 
%      probably because they came in that way from previous calculations
R1_Q_dallen = axisAngle2quaternion ( R1_aa_dallen )
% R1_Q_rotmat2quat = rotmat2quat( R1 ) % NOTE: matches this online algorithm if I use correct AA angle from vrrotmat2vec()
% R1_Q_dcm2quat = dcm2quat(R1) % check against this MATLAB function value
%      -- not available with student license? - https://www.mathworks.com/programs/trials/trial_request.html?prodcode=AT&eventid=572392830&s_iid=main_trial_AT_cta2
R2_Q_dallen = axisAngle2quaternion ( R2_aa_dallen )
% R2_Q_rotmat2quat = rotmat2quat( R2 ) % NOTE: matches this online algorithm if I use correct AA angle from vrrotmat2vec()
% R2_Q_dcm2quat = dcm2quat( R2 ) % check against this MATLAB function value 
%      -- not available with student license? - https://www.mathworks.com/programs/trials/trial_request.html?prodcode=AT&eventid=572392830&s_iid=main_trial_AT_cta2

%% slerp between the two quaternions
Qslerped_dallen = slerpQuaternions ( R1_Q_dallen, R2_Q_dallen, w )
% Qslerped_Dayot = slerpDayot ( R1_Q_dallen, R2_Q_dallen, w ) % NOTE: matches this online algorithm value

%% convert quaternion to axisAngle
AAslerped_dallen = quaternion2axisAngle ( Qslerped_dallen )
%      TODO: is there a matlab function for the above to check this? I currently
%      have no way to check if this is calculating properly

%% convert axisAngle to rotation matrix
%      NOTES: doesn't quite match the online algorythm yet, but that's
%      likely because my conversion from Quat to axis angle is off since my
%      calculations for this function match vrrotvec2mat() at this point.
Rslerped_dallen = axisAngle2rotmat ( AAslerped_dallen )
% Rslerped_vrrotvec2mat = vrrotvec2mat( AAslerped_dallen ) % NOTE: matches this MATLAB function value
Rslerped_quat2dc = quat2dcmTursa ( Qslerped_dallen ) % check against this online algorithm value
% R1_Q_quat2dcm = quat2dcm( Qslerped ) % check against this MATLAB function value
%      -- not available with student license? - https://www.mathworks.com/programs/trials/trial_request.html?prodcode=AT&eventid=572392830&s_iid=main_trial_AT_cta2
end

%% my own attempts at conversion functions

function [ AA ] = rotmat2axisAngle (R)
%ROTMAT2AXISANGLE Summary
%   Detailed explanation goes here

% get eigen vectors and values
[ V, D ] = eig(R);

% set axis to eigenVector column with corresponding value of 1
[ row, col ] = find(abs(1-D) < 0.0001); % TODO: less wasteful ways to do this? (lamda moded by size(D,2))?
rn = V(:, col);

% find theta 
%      TODO: angle number still incorrect
cosTheta = (trace(R)-1) / 2;
sincTheta = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]' ./ (2*rn);
sinTheta = sincTheta * norm(rn);

theta = theta(1); % NOTE: atan2 above returns a 3x1 matrix which ok since they're all the same, but we need a single number
theta = atan2( sinTheta, cosTheta ); 

AA = [rn; theta]; % TODO: signs are reversed (should technically be ok as long as all the numbers are correct)
end

function [ Q ] = axisAngle2quaternion ( AA )
%AXISANGLE2QUATERNION Summary
%   Detailed explanation goes here

Q = zeros([1, 4]);
theta = AA(4);
Q(1:3) = AA(1:3)' * sin(theta / 2);
Q(4) = cos(theta / 2);
end

function [ Qslerped ] = slerpQuaternions ( Q1, Q2, W )
%SLERP Summary
%   Detailed explanation goes here

theta = acos(dot(Q1, Q2));
Qslerped = (Q1 * sin((1-W) * theta) / sin(theta)) + (Q2 * sin(W * theta) / sin(theta));
end

function [ AA ] = quaternion2axisAngle ( Q )
%QUATERNION2AXISANGLE Summary
%   Detailed explanation goes here
theta = 2 * acos(Q(4));
rn = zeros([3, 1]);
rn(1:3) = Q(1:3)' ./ sin(theta / 2);

AA = [ rn; theta ];
end

function [ R ] = axisAngle2rotmat ( AA )
%AXISANGLE2ROTMAT Summary
%   Detailed explanation goes here
theta = AA(4);
rn = [ AA(1); AA(2); AA(3) ];
r = rn * theta;
rMag = norm(r); % TODO: doesn't this just equal theta?
rx = [                          ...
    0       -r(3)	r(2)    ;   ...
    r(3)    0       -r(1)   ;   ...
    -r(2)   r(1)    0           ...
    ];

R = cos(rMag) * eye(3)                  ...
    + sin(rMag) / rMag .* rx            ...
    + ((1-cos(rMag)) / rMag^2) * (r * r') ;
end

%% preexisting examples

% source - http://smallsats.org/2012/12/09/euler-rotation-example-rotation-matrix-quaternion-euler-axis-and-principal-angle/
function [ Q ] = rotmat2quat (R)
%ROTMAT2QUAT Summary
%
q4 = 0.5*(1 + R(1,1)+ R(2,2) + R(3,3))^0.5;
q1 = (R(2,3) - R(3,2))/(4*q4);
q2 = (R(3,1) - R(1,3))/(4*q4);
q3 = (R(1,2) - R(2,1))/(4*q4);

Q = [q1 q2 q3 q4];
end

% source - http://www.mathworks.com/matlabcentral/fileexchange/11827-slerp/content/slerp.m
%%%%%%%%%%%%%%%%%
%%%%     SLERP     %%%%%%
%%%%%%%%%%%%%%%%%

%        Sagi Dalyot %

%       This routine aims for calculating a unit quaternion,  describing a rotation matrix,
%       which lies between two known unit quaternions - q1 and q2,
%       using a spherical linear interpolation - Slerp.
%       Slerp follow the shortest great arc on a unit sphere,
%       hence, the shortest possible interpolation path.
%       Consequently, Slerp has constant angular velocity, 
%       so it is the optimal interpolation curve between two rotations.
%       (first published by Sheomake K., 1985 - Animating Rotation with Quaternion Curves)

%       end of file ->  explnation of rotation matrix and quaternions

%       in general:
%       slerp(q1, q2, t) = q1*(sin(1-t)*teta)/sin(t) + q2*(sin(t*teta))/sin(teta)
%       where teta is the angle between the two unit quaternions,
%       and t is between [0,1]

%       two border cases will be delt:
%       1: where q1 = q2 (or close by eps)
%       2: where q1 = -q2 (angle between unit quaternions is 180 degrees).
%       in general, if q1=q2 then Slerp(q; q; t) == q

function [qm] = slerpDayot (qi, qn, t, eps)

%       where qi=[w1 x1 y1 z1] - start unit quaternions
%                      qn=[w2 x2 y2 z2] - end unit quaternions
%                      t=[0 to 1]
%                      eps=threshold value

% added by dallen: set default epsilon if none was passed
if (nargin < 4)
    eps = 0.01;
end

if t==0 % saving calculation time -> where qm=qi
    qm=qi;
    
elseif t==1 % saving calculation time -> where qm=qn
    qm=qn;
    
else

    C=dot(qi,qn);                  % Calculating the angle beteen the unit quaternions by dot product

    teta=acos(C);

        if (1 - C) <= eps % if angle teta is close by epsilon to 0 degrees -> calculate by linear interpolation
            qm=qi*(1-t)+qn*t; % avoiding divisions by number close to 0

        elseif (1 + C) <= eps % when teta is close by epsilon to 180 degrees the result is undefined -> no shortest direction to rotate
            q2(1) = qi(4); q2(2) = -qi(3); q2(3)= qi(2); q2(4) = -qi(1); % rotating one of the unit quaternions by 90 degrees -> q2
            qm=qi*(sin((1-t)*(pi/2)))+q2*sin(t*(pi/2));

        else
            qm=qi*(sin((1-t)*teta))/sin(teta)+qn*sin(t*teta)/sin(teta);
        end
end
end
% eof
%  q = [w, (x, y, z)]    quaternion definition
% 
%  where
%  R = [1-2*y^2-2*z^2   2*x*y-2*w*z     2*x*z+2*w*y         R is function of 3  euler rotation angles
%         2*x*y+2*w*z    1-2*x^2-2*z^2   2*y*z-2*w*x
%         2*x*z-2*w*y     2*y*z+2*w*x    1-2*x^2-2*y^2]
%  => R = [M00 M01 M02
%                M10 M11 M12
%                M20 M21 M22]
%  => trace(R) = 4 - 4*(x^2+y^2+z^2), and x^2+y^2+z^2+w^2==1
%  => w=(trace)^.5
%  => x=(M21-M12)/4*w
%  => y=(M02-M21)/4*w
%  => x=(M10-M01)/4*w
%  => q = [w, (x, y, z)]


% source - http://www.mathworks.com/matlabcentral/newsreader/view_thread/160945
function [ dc ] = quat2dcmTursa(q)
% quat2dc quaternion direction cosine matrix angle axis
%*******************************************************************
%
% quat2dc calculates the dirction cosine matrix corresponding to a
% quaternion. Assumes input quaternion is normalized.
%
% Input: q = quaternion, q(1) = scalar, q(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
%
% Output: dc = 3x3 direction cosine matrix
%
% Programmer: James Tursa
%
%*******************************************************************

q11 = q(1)^2;
q12 = q(1)*q(2);
q13 = q(1)*q(3);
q14 = q(1)*q(4);
q22 = q(2)^2;
q23 = q(2)*q(3);
q24 = q(2)*q(4);
q33 = q(3)^2;
q34 = q(3)*q(4);
q44 = q(4)^2;

dc = zeros(3);
dc(1,1) = q11 + q22 - q33 - q44;
dc(2,1) = 2 * (q23 - q14);
dc(3,1) = 2 * (q24 + q13);
dc(1,2) = 2 * (q23 + q14);
dc(2,2) = q11 - q22 + q33 - q44;
dc(3,2) = 2 * (q34 - q12);
dc(1,3) = 2 * (q24 - q13);
dc(2,3) = 2 * (q34 + q12);
dc(3,3) = q11 - q22 - q33 + q44;

return;
end