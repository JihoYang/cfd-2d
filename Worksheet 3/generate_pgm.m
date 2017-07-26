clear all; clc; close all;

[folder, ~, ~] = fileparts(which('generate_pgm'))

C_B = 0;
C_F = 1;

% Box

imax = 100;
jmax = 20;

domain = ones(jmax, imax);



% % Karman Vortex Street
% domain(9:10, 12) = C_B;
% domain(9:11, 11) = C_B;
% domain(10:12, 10) = C_B;
% domain(11:12, 9) = C_B;
% imshow(domain);
% imwrite(domain, [folder, '/karman_vortex.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% Square object
% domain (10:15, 10:15) = C_B;
% imshow(domain);
% imwrite(domain, [folder, '/square.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% % Flow over step
% domain(jmax/2+1:end, 1:jmax/2) = C_B;
% imshow(domain);
% imwrite(domain, [folder, '/step.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);
% imwrite(domain, [folder, '/step.pgm'], 'ASCII', 'maxvalue', 1);

% % Plane shear flow
% imwrite(domain, [folder, '/plane_shear_flow.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% circle

%domain(5:7, 10) = C_B;
%domain(6, 9:11) = C_B;
%imshow(domain);
%imwrite(domain, [folder, '/circle.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);


