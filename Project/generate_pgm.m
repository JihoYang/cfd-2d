%-----------------------------------%
%                                   %
%   CFD Lab 2016: PGM generator     %
%                                   %
%   Group 8: Mohammed Asif Chand    %
%            Jiho Yang              %
%                                   %
%   Final Update: 12/07/2016        %
%                                   %
%-----------------------------------%

clear all; clc; close all;

[folder, ~, ~] = fileparts(which('generate_pgm'))

C_B = 0;
C_F = 1;

%Domain
imax = 100;
jmax = 60;

domain = ones(jmax, imax);


%-----------------------------------Geometries-------------------------------------

%----------------------------Multibody squares- layout 5-------------------

%First line
%Body 0
%domain(26:28, 20:22) = C_B;

%Second line
%Body 1
%domain(38:40, 40:42) = C_B;
%domain(13:15, 40:42) = C_B;
%domain(26:28, 40:42) = C_B;
%Body 2
%domain(26:28, 60:62) = C_B;

%Third line
%Body 3
%imshow(domain);
%imwrite(domain, [folder, '/layout5.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);



%---------------------------Multibody squares - layout 4---------------------------
%First line
%Body 0
%domain(26:28, 20:22) = C_B;

%Second line
%Body 1
%domain(38:40, 40:42) = C_B;
%domain(13:15, 40:42) = C_B;

%Body 2
%domain(26:28, 60:62) = C_B;

%Third line
%Body 3
%imshow(domain);
%imwrite(domain, [folder, '/layout4.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);
%---------------------------Multibody squares - layout 3---------------------------
%First line
%Body 0
%domain(30:32, 20:22) = C_B;

%Second line
%Body 1
%domain(40:42, 40:42) = C_B;
%domain(20:22, 40:42) = C_B;

%Body 2
%domain(30:32, 60:62) = C_B;

%Third line
%Body 3
%imshow(domain);
%imwrite(domain, [folder, '/layout3.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);


%---------------------------Multibody squares - layout 2---------------------------
%First line (only one line here)
%Body 0
%domain(50:52, 20:22) = C_B;

%Body 1
%domain(40:42, 20:22) = C_B;

%Body 2
%domain(30:32, 20:22) = C_B;

%Body 3
%domain(20:22, 20:22) = C_B;

%Body 4
%domain(10:12, 20:22) = C_B;

%imshow(domain);
%imwrite(domain, [folder, '/layout2.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

%--------------------------Multibody squares - layout 1---------------------------
%First line
%Body 0
%domain(50:52, 20:22) = C_B;

%Body 1
%domain(40:42, 20:22) = C_B;

%Body 2
%domain(30:32, 20:22) = C_B;

%Body 3
%domain(20:22, 20:22) = C_B;

%Body 4
%domain(10:12, 20:22) = C_B;

%Second line
%Body 5
%domain(45:47, 40:42) = C_B;

%Body 6
%domain(35:37, 40:42) = C_B;

%Body 7
%domain(25:27, 40:42) = C_B;

%Body 8
%domain(15:17, 40:42) = C_B;

%Third line
%Body 9
%domain(50:52, 60:62) = C_B;

%Body 10
%domain(40:42, 60:62) = C_B;

%Body 11
%domain(30:32, 60:62) = C_B;

%Body 12
%domain(20:22, 60:62) = C_B;

%Body 13
%domain(10:12, 60:62) = C_B;

%imshow(domain);
%imwrite(domain, [folder, '/layout1.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);


%---------------------------------Old Geometries-------------------------------

%Multibody circles - layout 1
% %Body 0
%domain(50, 20:21) = C_B;
%domain(51, 19:22) = C_B;
%domain(52, 19:22) = C_B;
%domain(53, 20:21) = C_B;
% 
% %Body 1
%domain(40, 20:21) = C_B;
%domain(41, 19:22) = C_B;
%domain(42, 19:22) = C_B;
%domain(43, 20:21) = C_B;
% 
% %Body 2
%domain(30, 20:21) = C_B;
%domain(31, 19:22) = C_B;
%domain(32, 19:22) = C_B;
%domain(33, 20:21) = C_B;
% 
% %Body 3
%domain(20, 20:21) = C_B;
%domain(21, 19:22) = C_B;
%domain(22, 19:22) = C_B;
%domain(23, 20:21) = C_B;
% 
% %Body 4
%domain(10, 20:21) = C_B;
%domain(11, 19:22) = C_B;
%domain(12, 19:22) = C_B;
%domain(13, 20:21) = C_B;
% 
% %Body 5
%domain(40, 50:51) = C_B;
%domain(41, 49:52) = C_B;
%domain(42, 49:52) = C_B;
%domain(43, 50:51) = C_B;
% 
% %Body 6
%domain(30, 50:51) = C_B;
%domain(31, 49:52) = C_B;
%domain(32, 49:52) = C_B;
%domain(33, 50:51) = C_B;
% 
% %Body 7
%domain(20, 50:51) = C_B;
%domain(21, 49:52) = C_B;
%domain(22, 49:52) = C_B;
%domain(23, 50:51) = C_B;
% 
%imshow(domain);
%imwrite(domain, [folder, '/multibody.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

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



