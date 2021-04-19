%Number of points required for a central difference of order d and truncation error O
%        p(d,O) = 2 * ⌊ (d+1) / 2 ⌋ - 1 + O
%
%Building a stencil:
%    The last term is constructed via truncation error order two elements,
%    the next previous one via error order four etc. until the second term
%    is reached. Afterwards all terms from 2 to n are subtracted from t1,
%    which is a padded version of the 3x3 stencil given below.
%
%                    1   0   1
%     3x3 stencil:   0  -4   0
%                    1   0   0
%
%    General construction of a stencil (up to the 13x13 case):
%        A13 = t1 - t2/6 - t3/180 - t4/10080 - t5/907200 - t6/119750400 - ...
%
%        t2 = dy4  +  6*dy2 *dx2 +         dx4
%        t3 = dy6  + 15*dy4 *dx2 +  15*dy2*dx4 +         dy6
%        t4 = dy8  + 28*dy6 *dx2 +  70*dy4*dx4 +  28*dy2*dx6 +         dx8
%        t5 = dy10 + 45*dy8 *dx2 + 210*dy6*dx4 + 210*dy4*dx6 +  45*dy2*dx8 +        dx10
%        t6 = dy12 + 66*dy10*dx2 + 495*dy8*dx4 + 924*dy6*dx6 + 495*dy4*dx8 + 66*dy2*dx10 + dx12

% Finite Difference Coefficients: {{{

% second derivative
d2o2  = [0, 0,    0,    0,     0,     1,     -2,     1,     0,    0,    0, 0, 0];
d2o4  = [0, 0,    0,    0,    -1,    16,    -30,    16,    -1,    0,    0, 0, 0]/12;
d2o6  = [0, 0,    0,    2,   -27,   270,   -490,   270,   -27,    2,    0, 0, 0]/180;
d2o8  = [0, 0,   -9,  128, -1008,  8064, -14350,  8064, -1008,  128,   -9, 0, 0]/5040;
d2o10 = [0, 8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8, 0]/25200;

dx2o2  = d2o2;  dy2o2  = d2o2';
dx2o4  = d2o4;  dy2o4  = d2o4';
dx2o6  = d2o6;  dy2o6  = d2o6';
dx2o8  = d2o8;  dy2o8  = d2o8';
dx2o10 = d2o10; dy2o10 = d2o10';

% fourth derivative
d4o2  = [  0,     0,     0,       0,       1,       -4,       6,       -4,       1,       0,     0,     0,   0];
d4o4  = [  0,     0,     0,      -1,      12,      -39,      56,      -39,      12,      -1,     0,     0,   0]/6;
d4o6  = [  0,     0,     7,     -96,     676,    -1952,    2730,    -1952,     676,     -96,     7,     0,   0]/240;
d4o8  = [  0,   -82,  1261,   -9738,   52428,  -140196,  192654,  -140196,   52428,   -9738,  1261,   -82,   0]/15120;
d4o10 = [479, -8208, 69444, -397520, 1809945, -4585248, 6222216, -4585248, 1809945, -397520, 69444, -8208, 479]/453600;

dx4o2  = d4o2;  dy4o2  = d4o2';
dx4o4  = d4o4;  dy4o4  = d4o4';
dx4o6  = d4o6;  dy4o6  = d4o6';
dx4o8  = d4o8;  dy4o8  = d4o8';
dx4o10 = d4o10; dy4o10 = d4o10';

% sixth derivative
d6o2  = [   0,     0,      0,      1,       -6,      15,      -20,      15,       -6,      1,      0,     0,    0];
d6o4  = [   0,     0,     -1,     12,      -52,     116,     -150,     116,      -52,     12,     -1,     0,    0]/4;
d6o6  = [   0,    13,   -190,   1305,    -4680,    9690,   -12276,    9690,    -4680,   1305,   -190,    13,    0]/240;
d6o8  = [-695, 11616, -93750, 481760, -1523385, 2992320, -3735732, 2992320, -1523385, 481760, -93750, 11616, -695]/60475;

dx6o2 = d6o2; dy6o2 = d6o2';
dx6o4 = d6o4; dy6o4 = d6o4';
dx6o6 = d6o6; dy6o6 = d6o6';
dx6o8 = d6o8; dy6o8 = d6o8';

% eighth derivative
d8o2  = [ 0,    0,    1,     -8,    28,    -56,    70,    -56,    28,     -8,    1,    0,  0];
d8o4  = [ 0,   -1,   13,    -69,   204,   -378,   462,   -378,   204,    -69,   13,   -1,  0]/3;
d8o6  = [31, -492, 3606, -15100, 39825, -69912, 84084, -69912, 39825, -15100, 3606, -492, 31]/360;

dx8o2 = d8o2; dy8o2 = d8o2';
dx8o4 = d8o4; dy8o4 = d8o4';
dx8o6 = d8o6; dy8o6 = d8o6';

% tenth derivative
d10o2 = [ 0,  1,  -10,   45, -120,  210,  -252,  210,  -120,   45,  -10,  1,  0];
d10o4 = [-5, 72, -450, 1640, 3915, 6480, -7644, 6480, -3915, 1640, -450, 72, -5]/12;

dx10o2 = d10o2; dy10o2 = d10o2';
dx10o4 = d10o4; dy10o4 = d10o4';

d12o2 = [1, -12, 66, -220, 495, -792, 924, -792, 495, -220, 66, -12, 1];

dx12o2 = d12o2; dy12o2 = d12o2';
%}}}

% All coefficient arrays should have the same size
arraysize  = 13;
centralidx = floor(arraysize/2) + 1;

% prepare the 3x3 stencil, which is also the term t1 for all higher order stencils
t1 = zeros(arraysize,arraysize);
t1( (arraysize-3)/2 + 1 : (arraysize-3)/2 + 3 , (arraysize-3)/2 + 1 : (arraysize-3)/2 + 3 ) = [1 0 1; 0 -4 0; 1 0 1];

% function to truncate stencil to its final size
function A = TruncateStencil(A, stencilsize)
    arraysize = rows(A);
    A = A( (arraysize-stencilsize)/2 + 1 : (arraysize-stencilsize)/2 + stencilsize , (arraysize-stencilsize)/2 + 1 : (arraysize-stencilsize)/2 + stencilsize );
endfunction

% function to write the stencil to a file
function writestencil(A, multiplier, filename)

    fd = fopen(filename, "w");

    for i = 1:prod(size(A))
        fprintf(fd, "stencil[%3d] = %f/%f;\n", i-1, A(i)*multiplier, multiplier);
    end

    fclose(fd);

endfunction

% Stencil generation
%--------------------------------------------------------------------------------
% 3x3 stencil {{{
%--------------------------------------------------
A3 = TruncateStencil(t1, 3);
%}}}
% 5x5 stencil {{{
%--------------------------------------------------
t2 = 6*dy2o2*dx2o2;         % 6/1
t2(centralidx,:) += dx4o2;  % 1/1
t2(:,centralidx) += dy4o2;  % 1/1
% lcm denominator of t2/6:    1/1

A5 = t1 - t2/6;
A5 = TruncateStencil(A5, 5);
%}}}
% 7x7 stencil {{{
%--------------------------------------------------
t2 = 6*dy2o4*dx2o4;         % 6/12/12
t2(centralidx,:) += dx4o4;  % 1/6
t2(:,centralidx) += dy4o4;  % 1/6
% lcm denominator of t2/6:    1/144

t3  = 15*dy2o2*dx4o2;       % 15/1
t3 += 15*dy4o2*dx2o2;       % 15/1
t3(centralidx,:) += dx6o2;  % /1
t3(:,centralidx) += dy6o2;  % /1
% lcm denominator of t3/180: 1/180

% lcm(144,180) = 720
A7 = t1 - t2/6 - t3/180;
A7 = TruncateStencil(A7, 7);
%}}}
% 9x9 stencil {{{
%--------------------------------------------------
t2 = 6*dy2o6*dx2o6;         % 6/180/180
t2(centralidx,:) += dx4o6;  % 1/240
t2(:,centralidx) += dy4o6;  % 1/240
% lcm denominator of t2/6:    1/32400

t3  = 15*dy2o4*dx4o4;       % 15/12/6
t3 += 15*dy4o4*dx2o4;       % 15/6/12
t3(centralidx,:) += dx6o4;  % 1/4
t3(:,centralidx) += dy6o4;  % 1/4
% lcm denominator of t3/180:  1/4320

t4  = 28*dy2o2*dx6o2;       % 28/1
t4 += 70*dy4o2*dx4o2;       % 70/1
t4 += 28*dy6o2*dx2o2;       % 28/1
t4(centralidx,:) += dx8o2;  %  1/1
t4(:,centralidx) += dy8o2;  %  1/1
% lcm denominator of t3/10080: 1/10080

% lcm(32400,4320,10080) = 453600
A9 = t1 - t2/6 - t3/180 - t4/10080;
A9 = TruncateStencil(A9, 9);
%}}}
% 11x11 stencil {{{
%--------------------------------------------------
t2 = 6*dy2o8*dx2o8;         % 6/5040/5040
t2(centralidx,:) += dx4o8;  % 1/15120
t2(:,centralidx) += dy4o8;  % 1/15120
% lcm denominator of t2/6:    1/25401600

t3  = 15*dy2o6*dx4o6;       % 15/180/240
t3 += 15*dy4o6*dx2o6;       % 15/240/180
t3(centralidx,:) += dx6o6;  %  1/240
t3(:,centralidx) += dy6o6;  %  1/240
% lcm denominator of t3/180:   1/518400

t4  = 28*dy2o4*dx6o4;       % 28/12/4
t4 += 70*dy4o4*dx4o4;       % 70/6/6
t4 += 28*dy6o4*dx2o4;       % 28/4/12
t4(centralidx,:) += dx8o4;  %  1/3
t4(:,centralidx) += dy8o4;  %  1/3
% lcm denominator of t4/10080: 1/362880

t5  =  45*dy2o2*dx8o2;      %  45/1
t5 += 210*dy4o2*dx6o2;      % 210/1
t5 += 210*dy6o2*dx4o2;      % 210/1
t5 +=  45*dy8o2*dx2o2;      %  45/1
t5(centralidx,:) += dx10o2; %   1/1
t5(:,centralidx) += dy10o2; %   1/1
% lcm denominator of t5/907200: 1/907200

% lcm(25401600,518400,362880,907200) = 25401600
A11 = t1 - t2/6 - t3/180 - t4/10080 - t5/907200;
A11 = TruncateStencil(A11, 11);
%}}}
% 13x13 stencil {{{
%--------------------------------------------------
t2 = 6*dy2o10*dx2o10;       % 6/25200
t2(centralidx,:) += dx4o10; % 1/453600
t2(:,centralidx) += dy4o10; % 1/453600
% lcm denominator of t2/6:    1/2721600

t3  = 15*dy2o8*dx4o8;       % 15/5040/15120
t3 += 15*dy4o8*dx2o8;       % 15/15120/5040
t3(centralidx,:) += dx6o8;  %  1/60475
t3(:,centralidx) += dy6o8;  %  1/60475
% lcm denominator of t3/180:   1/11060364672000

t4  = 28*dy2o6*dx6o6;       % 28/180/240
t4 += 70*dy4o6*dx4o6;       % 70/240/240
t4 += 28*dy6o6*dx2o6;       % 28/240/180
t4(centralidx,:) += dx8o6;  %  1/360
t4(:,centralidx) += dy8o6;  %  1/360
% lcm denominator of t4/10080: 1/870912000

t5  =  45*dy2o4*dx8o4;      %  45/12/3
t5 += 210*dy4o4*dx6o4;      % 210/6/4
t5 += 210*dy6o4*dx4o4;      % 210/4/6
t5 +=  45*dy8o4*dx2o4;      %  45/3/12
t5(centralidx,:) += dx10o4; %   1/12
t5(:,centralidx) += dy10o4; %   1/12
% lcm denominator of t5/907200: 1/10886400

t6  =  66*dy2o2*dx10o2;     %  66/1
t6 += 495*dy4o2*dx8o2;      % 495/1
t6 += 924*dy6o2*dx6o2;      % 924/1
t6 += 495*dy8o2*dx4o2;      % 495/1
t6 +=  66*dy10o2*dx2o2;     %  66/1
t6(centralidx,:) += dx12o2; %   1/1
t6(:,centralidx) += dy12o2; %   1/1
% lcm denominator of t6/119750400:  1/119750400

% lcm(2721600,11060364672000,870912000,10886400,119750400) = 486656045568000
A13 = t1 - t2/6 - t3/180 - t4/10080 - t5/907200 - t6/119750400;
A13 = TruncateStencil(A13, 13);
%}}}
