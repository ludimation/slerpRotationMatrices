function [ Rslerped_dallen ] = slerpRotationMatrices( R1, R2, w )
%SLERPROTATIONMATRICES Summary
%   Detailed explanation goes here
%      TODO: might find some more ideas here ? http://www.mathworks.com/matlabcentral/newsreader/view_thread/301137
%      -- some other methods not available with student license? - https://www.mathworks.com/programs/trials/trial_request.html?prodcode=AT&eventid=572392830&s_iid=main_trial_AT_cta2
%         R2_Q_dcm2quat = dcm2quat( R2 ); % check against this MATLAB function value 
%         R1_Q_quat2dcm = quat2dcm( Qslerped ); % check against this MATLAB function value

%% convert rotation matrices to axisAngles
R1_aa_dallen = rotMat2axisAngle ( R1 )
R2_aa_dallen = rotMat2axisAngle ( R2 )
% % NOTE: the values returned by functions above matche values ouput by
%      MATLAB functions below
% R1_aa_vrrotmat2vec = vrrotmat2vec( R1 )
% R2_aa_vrrotmat2vec = vrrotmat2vec( R2 )


%% convert axisAngles to quaternions
R1_Q_dallen = axisAngle2quaternion ( R1_aa_dallen )
R2_Q_dallen = axisAngle2quaternion ( R2_aa_dallen )
% % NOTE: values from above seem correct since they match values from an
%      online algorithm (below), but signs are reversed. is this where my
%      problems are coming from?
% R1_Q_rotmat2quat = rotmat2quat( R1 ) % NOTE: matches this online algorithm
% R2_Q_rotmat2quat = rotmat2quat( R2 ) % NOTE: matches this online algorithm

%% slerp between the two quaternions
Qslerped_dallen = slerpQuaternions ( R1_Q_dallen, R2_Q_dallen, w )
% % NOTE: value returned above matches value of online algorithm below
% Qslerped_Dayot = slerpDayot ( R1_Q_dallen, R2_Q_dallen, w )

%% convert quaternion to axisAngle
AAslerped_dallen = quaternion2axisAngle ( Qslerped_dallen )
% % NOTE: value returned above matches the value returned by the inverse 
%      function below
% Qslerped_dallen_axisAngle2quaternion = axisAngle2quaternion( AAslerped_dallen )

%% convert axisAngle to rotation matrix
Rslerped_dallen = axisAngle2rotMat ( AAslerped_dallen )
% % NOTE: value does not match the online algorythm Rslerped_quat2dc(), or
%      the MATLAB function vrrotvec2mat(), nor the inverse functions
%      rotMat2axisAngle() and MATLAB's vrrotmat2vec()
Rslerped_vrrotvec2mat = vrrotvec2mat( AAslerped_dallen )
Rslerped_quat2dc = quat2dcmTursa ( [ Qslerped_dallen(4), Qslerped_dallen(1:3) ] )    % assumes scalar is in the first position of the matrix
% % test inverse funtions -- values below should equal AAslerped_dallen (not working yet :()
% AAslerped_dallen_dallenRotMat2axisAngle = rotMat2axisAngle ( Rslerped_dallen )
% AAslerped_dallen_vrrotmat2vec = vrrotmat2vec ( Rslerped_dallen )
end

%% my own attempts at conversion functions

function [ AA ] = rotMat2axisAngle ( R )
%ROTMAT2AXISANGLE Summary
%   Detailed explanation goes here

% get eigen vectors and values
[ V, D ] = eig( R );

% set axis to eigenVector column with corresponding value of 1
[ row, col ] = find( abs( 1 - D ) < 0.0001 ); % TODO: less wasteful ways to do this? (lamda moded by size(D,2))?
rn = V( :, col );

% find theta 
cosTheta = ( trace( R ) - 1 ) / 2;
sincTheta = [ R( 3, 2 ) - R( 2, 3 ), R( 1, 3 ) - R( 3, 1 ), R( 2, 1 )-R( 1, 2 ) ]' ./ ( 2 * rn );
sinTheta = sincTheta * norm( rn );

theta = atan2( sinTheta, cosTheta ); 
theta = theta( 1 ); % NOTE: atan2 above returns a 3x1 matrix which is ok since all the values are the same, but we need a single number

AA = [ rn; theta ]; % TODO: signs are reversed (should technically be ok since it is simmetrical to negated axis angle... as long as all the numbers are correct)
end

function [ Q ] = axisAngle2quaternion ( AA )
%AXISANGLE2QUATERNION Summary
%   assumes AA = axis angle representation with AA(4) = theta, AA(1:3) = vector
%   returns Q with Q(4) = scalar, Q(1:3) = vector (is this a normalized quaterion?)

Q = zeros( [ 1, 4 ] );
theta = AA( 4 );
Q( 1:3 ) = AA( 1:3 )' * sin( theta / 2 );
Q( 4 ) = cos( theta / 2 );
end

function [ Qslerped ] = slerpQuaternions ( Q1, Q2, W )
%SLERP Summary
%   Detailed explanation goes here

theta = acos( dot( Q1, Q2 ) );
Qslerped = ( Q1 * sin( ( 1 - W ) * theta ) / sin( theta ) ) + ( Q2 * sin( W * theta ) / sin( theta ) );
end

function [ AA ] = quaternion2axisAngle ( Q )
%QUATERNION2AXISANGLE Summary
%   assumes Q = Quaterion representation with Q(4) = scalar, Q(1:3) = vector (is this a normalized quaterion?)
%   returns AA = axis angle representation with AA(4) = theta, AA(1:3) = vector
theta = 2 * acos( Q( 4 ) );
rn = zeros( [ 3, 1 ] );
rn( 1:3 ) = Q( 1:3 )' ./ sin( theta / 2 );

AA = [ rn; theta ];
end

function [ R ] = axisAngle2rotMat ( AA )
%AXISANGLE2ROTMAT Summary
%   Detailed explanation goes here
theta = AA( 4 );
rn = AA ( 1:3 );
r = rn * theta;
rMag = abs( theta );
% create rx based on the scaled r (magnitude of theta)?
rx = [                        ...
    0       -r(3)	r(2)	; ...
    r(3)	0       -r(1)	; ...
    -r(2)	r(1)	0         ...
    ];
% or should I create rx based on the normalized rn (magnitude of 1)?
rxn = [                        ...
    0       -rn(3)	rn(2)   ; ...
    rn(3)   0       -rn(1)	; ...
    -rn(2)  rn(1)   0         ...
    ];

secondTermSinc = sinc(rMag) * rxn;
secondTermSinOverMag = (sin(pi * rMag) / (pi * rMag)) * rxn;

% TODO: this calculation seems to still be a bit off (diagonal values are
%      correct, but all upper right and lower left values are not. Does
%      that suggest there is something wrong with the rx matrix or the
%      second term in the function below?)
% R = cos(rMag) * eye(3)                              ...    
%     + sinc(rMag) * rx                               ... (sin(pi * rMag) / (pi * rMag)) * rx
%     + ( ( 1 - cos(rMag) ) / rMag^2 ) * ( r * r' )   ; 
% another version of Rodriguez formula from wikipedia -- http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula#Conversion_to_rotation_matrix
R = cos(rMag) * eye(3)           	...    
    + sin(rMag) * rxn               	...
    + ( 1 - cos(rMag) ) * ( rn * rn' )	; 
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