function Visual_flow(Bn_se, Bn_mo)

% Bn_se = 1;
% Bn_mo = 1;
N_mode = length(Bn_se(1,:)); % number of velocity modes

N = 400; % number of grids for computing fluid
zrange = 5; xrange = 5; % computation domain in (z,x) plane
z = - zrange : 2*zrange/N : zrange;
x = -xrange:2*xrange/N:xrange; 
x0 = x(N/2 + 1:end)'; % x >= 0 
[XX,ZZ] = meshgrid(x,z); % upper half plane 

% convert to spherical coordinates
r = sqrt(z.^2 + x0.^2);
sintheta = x0./r;
costheta = z./r;

% require space for velocity data
ur0 = zeros(length(x0), length(z)); 
ut0 = zeros(length(x0), length(z));

ure = ur0;  ute = ut0; 
urm = ur0;  utm = ut0; 

P_mu = costheta; P_mumu = costheta; % Legendre polynomial n=1
for n = 1 : N_mode
    [fre, fte] = vmode_sessile(r, n); % Sessile modes   
    [frm, ftm] = vmode_swim(r, n); % Motile modes 
    [Pn_new, Vn_new, P_mu, P_mumu] = PnVn_cart(costheta, n, P_mu, P_mumu);
    % flow around sessile cell
    ure = ure + Bn_se(n).*fre.*Pn_new; 
    ute = ute + Bn_se(n).*fte.*Vn_new; 
    % flow around motile cell
    urm = urm + Bn_mo(n).*frm.*Pn_new; 
    utm = utm + Bn_mo(n).*ftm.*Vn_new;
end 

% convert velocity into cartesian coordinates
[uze0, uxe0] = spherical_to_cartesian(ure, ute, costheta, sintheta);
[uzm0, uxm0] = spherical_to_cartesian(urm, utm, costheta, sintheta);
% lab frame velocity for motile
uzmlab0 = uzm0 + (2/3)*Bn_mo(1);
uxmlab0 = uxm0;

% recove full velocity field using phi-axis symmetry
[uze, uxe] = recover_full_u(uze0, uxe0, r);
[uzm, uxm] = recover_full_u(uzm0, uxm0, r);
[uzmlab, uxmlab] = recover_full_u(uzmlab0, uxmlab0, r);

%% visual fluid in streamslice
cob = [0, 0.4470, 0.7410]; % blue
theta = 0:pi/100:2*pi;
Circle1 = cos(theta);   Circle2 = sin(theta); 

bound = 5; % set plot boundary
figure('color','w')
subplot(2,3,1) % sessile
h1 = streamslice(XX, ZZ, uxe',uze', 0.5, 'linear');
set(h1, 'LineWidth',1,'color',cob)
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k')
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('sessile')

subplot(2,3,2) % motile 
h2 = streamslice(XX, ZZ, uxm', uzm', 0.5, 'linear');
set(h2, 'LineWidth',1,'color',cob)
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k')
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('motile (body frame)')

subplot(2,3,3) % motile 
h3 = streamslice(XX, ZZ, uxmlab', uzmlab', 0.5, 'linear');
set(h3, 'LineWidth',1,'color',cob)
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k')
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('motile (lab frame)')

%% velocity magnitude
cob = [0, 0.4470, 0.7410]; % blue
[grad,im]=colorGradient([1,1,1],cob,256);

uemag = sqrt(uze.^2+uxe.^2); ue_max = max(max(uemag));
ummag = sqrt(uzm.^2+uxm.^2); um_max = max(max(ummag));
umlabmag = sqrt(uzmlab.^2+uxmlab.^2); umlab_max = max(max(umlabmag));

subplot(2,3,4)
pcolor(XX,ZZ,uemag'); 
colormap(grad); caxis([0,ue_max]); shading interp; colorbar
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k') 
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('sessile')

subplot(2,3,5)
pcolor(XX,ZZ,ummag'); 
colormap(grad); caxis([0,um_max]); shading interp; colorbar
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k')
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('motile (body frame)')

subplot(2,3,6)
pcolor(XX,ZZ,umlabmag'); 
colormap(grad); caxis([0,umlab_max]); shading interp; colorbar
hold on; plot(Circle1,Circle2, 'linewidth',1, 'color', 'k')
axis equal; xlim([-bound,bound]); ylim([-bound,bound])
title('motile (lab frame)')


% Use recursion relation to calculate Legnedre polynomials
function [Pn_new, Vn_new, P_mu, P_mumu] = PnVn_cart(co, n, P_mu, P_mumu)
    if n == 1
        Pn_new = co;
        Vn_new = sqrt(1-co.^2);
        P_mumu = ones(length(co(:,1)),length(co(1,:)));
    else
        Pn_new = ( (2*n-1).*co.*P_mu - (n-1).*P_mumu )./n;
        dP_new = ( co.*Pn_new - P_mu ).*n./( co.^2-1 );
        Vn_new = sqrt(1-co.^2).*dP_new.*2./(n*(n+1));
        Vn_new(co == 1 ) = 0; Vn_new(co == -1 ) = 0;
        P_mumu = P_mu;
    end
    P_mu = Pn_new;
    Vn_new(:,1) = 0;    
    end

function [fr, ft] = vmode_sessile(r, n) % Sessile modes

    fr = (1./r.^(n+2) - 1./r.^(n));
    ft = (n./r.^(n+2) - (n-2)./r.^(n))./2;  

end

function [fr, ft] = vmode_swim(r, n) % Motile modes

    if n == 1
        fr = ( 1./r.^3 - 1 ).*2/3;
        ft = ( 1./r.^3 + 2 )./3;
    else
        fr = 1./r.^(n+2) - 1./r.^(n);
        ft = (n./r.^(n+2) - (n-2)./r.^(n))./2;  
    end

end

% convert flow field into cartesian coordinates
function [uz, ux] = spherical_to_cartesian(ur, ut, costheta, sintheta)

    uz = ur.*costheta - ut.*sintheta;  
    ux = ur.*sintheta + ut.*costheta;

end 


function [uz, ux] = recover_full_u(uz0, ux0, r)
    % remove calculation inside sphere
    uz0(r<1) = 0; ux0(r<1) = 0;
    % recover full field
    uz = [flip(uz0(2:end,:),1); uz0];
    ux = [-flip(ux0(2:end,:),1); ux0];
    
end


function [grad,im]=colorGradient(c1,c2,depth)
% COLORGRADIENT allows you to generate a gradient between 2 given colors,
% that can be used as colormap in your figures.
%
% USAGE:
%
% [grad,im]=getGradient(c1,c2,depth)
%
% INPUT:
% - c1: color vector given as Intensity or RGB color. Initial value.
% - c2: same as c1. This is the final value of the gradient.
% - depth: number of colors or elements of the gradient.
%
% OUTPUT:
% - grad: a matrix of depth*3 elements containing colormap (or gradient).
% - im: a depth*20*3 RGB image that can be used to display the result.
%
% EXAMPLES:
% grad=colorGradient([1 0 0],[0.5 0.8 1],128);
% surf(peaks)
% colormap(grad);
%
% --------------------
% [grad,im]=colorGradient([1 0 0],[0.5 0.8 1],128);
% image(im); %display an image with the color gradient.

% Copyright 2011. Jose Maria Garcia-Valdecasas Bernal
% v:1.0 22 May 2011. Initial release.

%Check input arguments.
%input arguments must be 2 or 3.
error(nargchk(2, 3, nargin));

%If c1 or c2 is not a valid RGB vector return an error.
if numel(c1)~=3
    error('color c1 is not a valir RGB vector');
end
if numel(c2)~=3
    error('color c2 is not a valir RGB vector');
end

if max(c1)>1&&max(c1)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c1 is not given as intensity values. Trying to convert');
    c1=c1./255;
elseif max(c1)>255||min(c1)<0
    error('C1 RGB values are not valid.')
end

if max(c2)>1&&max(c2)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c2 is not given as intensity values. Trying to convert');
    c2=c2./255;
elseif max(c2)>255||min(c2)<0
    error('C2 RGB values are not valid.')
end
%default depth is 64 colors. Just in case we did not define that argument.
if nargin < 3
    depth=64;
end

%determine increment step for each color channel.
dr=(c2(1)-c1(1))/(depth-1);
dg=(c2(2)-c1(2))/(depth-1);
db=(c2(3)-c1(3))/(depth-1);

%initialize gradient matrix.
grad=zeros(depth,3);
%initialize matrix for each color. Needed for the image. Size 20*depth.
r=zeros(20,depth);
g=zeros(20,depth);
b=zeros(20,depth);
%for each color step, increase/reduce the value of Intensity data.
for j=1:depth
    grad(j,1)=c1(1)+dr*(j-1);
    grad(j,2)=c1(2)+dg*(j-1);
    grad(j,3)=c1(3)+db*(j-1);
    r(:,j)=grad(j,1);
    g(:,j)=grad(j,2);
    b(:,j)=grad(j,3);
end

%merge R G B matrix and obtain our image.
im=cat(3,r,g,b);

end
end