clc
close all;
clear all;
pkg load image;
strp = './data/input/'; % % where the images are placed
strpy = './data/output/'; %ylocation
%strpl='/Users/debjanisingharoy/Desktop/Suman/Suman Research details/research proposal/yarin/Droplet_machine learning/Motion for Suman/analyzed files/location/'; %location
str = 'water, 250 ms, 0.148 kV '; %name of files without counter
extstr = 'jpg'; %extension of image file
k = 0; % % this is a counter but we don't need it at later time

%%**** define the scale in the system*******************************************
%d_scale = 1.0; % calibration ball diameter [mm]
%A= imread('/Users/debjanisingharoy/Desktop/Suman/Suman Research details/research proposal/yarin/Droplet_machine learning/Motion for Suman/Scale.jpg'); %%% Scale bar
%imshow(A);
%title('MAKE CALIBRATION !!');
%[x, y]  = ginput(1);
%v(1, :) = [x, y];
%[x, y]  = ginput(1);
%v(2, :) = [x, y];
%factor = d_scale/(v(2, 1) - v(1, 1));
factor = 0.0048742;
k1 = 0;
%%%%****************************************************************************

%% input parameters
%first= input("what is the first image in the series -->  "); %%add the first value
%last= input("what is the last image in the series -->  "); %%add the last value
first = 1;
last = 1100;
step = 5; % %provides step between analyzed images
fps = 1500; %frames per second
level = 2; % %use this to define loop
n = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = first:step:last % Main for loop
    %%%%% This loop tells you the image counter%%%%%%%%
    if (i < 10),
        str1 = '000';
        str2 = [str, num2str(str1), num2str(i)];

    elseif ((i >= 10) && (i < 100)),
        str1 = '00'; str2 = [str, num2str(str1), num2str(i)];

    elseif ((i >= 100) && (i < 1000)),
        str1 = '0'; str2 = [str, num2str(str1), num2str(i)];

    else
        str2 = [str, num2str(i)];

    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([strp, str2, '.', extstr]);
    A = imread([strp, str2, '.', extstr], extstr); % % read the image

    k = 0;
    n = n + 1; % %counter for sub-models
    %imshow(A)
    if level == 1
        [a1, a2, imc, siz] = imcrop(A); % %crop the image to reasonable size & find the co-ordinates
        siz = round(siz); % % Provides co-ordinates of the rectangle edges

        %%Defines all the co-ordinates************
        x1 = siz(1);
        y1 = siz(2);
        x2 = siz(1) + siz(3);
        y2 = siz(2) + siz(4);
        %level = 0;
    end;

    %%%%***************************************
    if level == 2
        x1 = 19;
        y1 = 419;
        x2 = 1900;
        y2 = 825;
        %level = 0;
    end;

    C = im2bw(A(y1:y2, x1:x2), graythresh(A(y1:y2, x1:x2))); % %change to b&W... 50 to avoid any unwanted specs
    A = double(C);
    %  imshow(A);
    [d1 d2] = size(A); % % computes the overall size
    T1(:, 1) = 0:(d2 - 1); % % X-co-ordinates
    %i %%%for debugging
    T1(:, 1) = T1(:, 1) .* factor;
    T1(:, 2) = d1;

    for i1 = 1:d2
        S1 = 1; % % In case of image, where the boundary is black start with 1 and for white start with 0. Opposite of intended boundary value needed
        S1 = find(A(:, i1) == 0); % % For black boundary, this is 0 and for white, this is 1
        S = min(S1); % % the first place where the boundary inversion happens

        if (S >= 1)
            T1(i1, 2) = S; % % store the value
        end;

    end;

    T1(:, 2) = d1 - T1(:, 2); % %final height at a certain x
    T1(:, 2) = T1(:, 2) .* factor;
    %%%%%%%%Level set and setup the origin at t = 0 and use that to study

    M = max(T1(:, 2)); % % %find the tip of the circular arc
    A = find(T1(:, 2) == M); % % % Find the places where all the plateau is reached
    tip = round((min(A) + max(A)) * 0.5); % % % Find the middlemost X co-ordinate value
    level = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%This section looks at the stretch of the droplet
    S4 = T1(:, 2);
    S5 = find(S4 > 0);
    lengthy = max(S5) - min(S5);
    length(n) = lengthy * factor;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tiporigin(n) = tip;
    timestamp(n) = i;

    for j1 = 1:50:d2
        k = k + 1;
        y(k) = T1(j1, 2);
        %X(k) = j1;
    end;

    y = y';

    for j1 = 1:k
        p1 = y(j1);
        yml(j1, n) = p1;
    end;

end;

savefile = strcat([[strpy, str], '_yml.dat']);
disp(savefile);
save(savefile, 'yml', '-ascii');
%M1= min(tiporigin);
%A1= find (tiporigin == M1);
%timeorigin = timestamp(A1);
%k=0;

time = time';
savefile = strcat([[strpy, str], '_time.dat']);
disp(savefile);
save(savefile, 'time', '-ascii');

%xposition = T1(:,1); %%% Define the x co-ordinates in a vector
M6 = max(length); % % Find where it is stretched the most
S12 = find (length == M6); % % Location of the highest value
xorigin = tiporigin (S12) * factor; % % Find the modified origin
timeorigin = timestamp(S12);

for j1 = 1:50:d2
    k1 = k1 + 1;
    xpos(k1) = T1(j1, 1);
    %X(k) = j1;
end;

for j1 = 1:k1
    x1 = xpos(j1);
    xcoordinate(j1) = x1 - xorigin;
end;


k = 0;
xcoordinate = xcoordinate';

for i = first:step:last
    k = k + 1;
    time(k) = i - timeorigin;

end

disp(time);
exit(1);

time = time';
savefile = strcat([[strpy, str], '_time.dat']);
disp(savefile);
save(savefile, 'time', '-ascii');
savefile = strcat([[strpy, str], '_xcoord.dat']);
disp(savefile);
save(savefile, 'xcoordinate', '-ascii');
