%% Main code for solving steadty sate advection-diffusion equation by using finite difference method for concentration field. 
% Fluid field is computed based on prescribed surface velocity on cell. 

clear; clc
% parameters ==============================================================
R = 20; % Computation region
Nr = 801; % grid number in r space
feeding_area = 0.5; % 50 percent feeding area
[theta, indcap] = theta_space(feeding_area);
[d1, d2, r, dzeta] = StretchMesh(R, Nr, 1 ,1); % stretch r space
Pe  = 100;

% Surface velocity and coefficients Bn ------------------------------------
N_mode = 100; % number of velocity modes, larger number of modes are needed for small fraction of slip velocity
[Bn_se, Bn_mo] = surface_velocity_expansion(feeding_area, -1, N_mode); % velocity coefficients in sessile and motile models

[urs, uts] = Flow_Sessile(r, theta, Bn_se, N_mode); % Flow field sessile
[urm, utm] = Flow_Motile(r, theta, Bn_mo, N_mode); % Flow field motile

%% plot fluid field arround model cell
Visual_flow(Bn_se, Bn_mo);
saveas(gcf, '../results/flow_fields.png')

%% Concentration field calculation =========================================
% pure diffusion, no flow
[c_0, I_0] = Advection_Diffusion_varycap_stretchmesh(r, d1, d2, dzeta, theta, urm, utm, 0, indcap);
% Advection-diffusion in motile cell
[c_mot, I_mot] = Advection_Diffusion_varycap_stretchmesh(r, d1, d2, dzeta, theta, urm, utm, Pe, indcap);
% Advection-diffusion in sessile cell
[c_ses, I_ses] = Advection_Diffusion_varycap_stretchmesh(r, d1, d2, dzeta, theta, urs, uts, Pe, indcap);

Sh_mot = I_mot/I_0;
Sh_ses = I_ses/I_0;

% Show concenration fields in body frame of reference
Visual_concentration(r, theta, c_ses,c_mot,Sh_ses,Sh_mot);
saveas(gcf, '../results/concentration_fields.png')

% Save to files, uncomment it if need to save data
% savename = ['Pe100_concentration_capture_50percent_cilia_50_percent.mat'];
% save(savename);

%% Functions
function [urm, utm] = Flow_Sessile(r, theta, Bn, Nv)
    [Pn, Vn] = PmVm_mu(cos(theta), Nv+1);
    
    rv = r;
    urm = zeros(length(r), length(theta));
    utm = zeros(length(r), length(theta));
    for n = 1:Nv
        urm = urm + (1./rv.^(n+2) - 1./rv.^(n)).*Bn(n)*Pn(n+1,:);
        utm = utm + 0.5.*(n./rv.^(n+2) - (n-2)./rv.^(n)).*Bn(n)*Vn(n+1,:);  
    end
end

function [urm, utm] = Flow_Motile(r, theta, Bn, Nv)
    Nr = length(r);
    
    [Pn, Vn] = PmVm_mu(cos(theta), Nv+1);
    
    rv = r;
    urm = (1./rv.^3 - ones(Nr,1)).*2/3.*Bn(1).*Pn(2,:);
    utm = (1./rv.^3 + 2*ones(Nr,1))./3.*Bn(1).*Vn(2,:);  
    
    if Nv > 1      
    for n = 2:Nv
        urm = urm + (1./rv.^(n+2) - 1./rv.^(n)).*Bn(n)*Pn(n+1,:);
        utm = utm + 0.5.*(n./rv.^(n+2) - (n-2)./rv.^(n)).*Bn(n)*Vn(n+1,:);  
    end
    end
end

% Legendre polynomails
function [P_mu, V_mu] = PmVm_mu(mu, LM)
    
    Ntheta = length(mu);
    P_mu = zeros(LM, Ntheta); 
    dP_mu = zeros(LM, Ntheta);
    V_mu = zeros(LM, Ntheta);
    
    P_mu(1,:) = 1; dP_mu(1,:) = 0; V_mu(1,:) = 0; % 0 mode
    P_mu(2,:) = mu; dP_mu(2,:) = 1; V_mu(2,:) = sqrt(1-mu.^2);% 1 mode
    
    for mm = 3:LM
        m = mm-1;
        P_mu(mm,:) = ((2*m-1)*mu.*P_mu(mm-1,:)-(m-1)*P_mu(mm-2,:))/m;
        dP_mu(mm,:) = (mu.*P_mu(mm,:) - P_mu(mm-1,:)).*m./(mu.^2-1);
        V_mu(mm,:) = sqrt(1-mu.^2).*dP_mu(mm,:)*2/(m*(m+1));
    end

end

function [theta, indcap] = theta_space(feeding_area)
% need adjust according to feeding area
    dtheta0 = pi/400;
    theta = 0:dtheta0:pi;
    newvec = (theta - acos(1-2.*feeding_area) ).^2; 
    % check if mu2 is accurate 
    % minnewvec = min(newvec) 
    indcap = find(newvec == min(newvec)); % index of finding mu2
    
end

function [d1, d2, r, dzeta] = StretchMesh(Rmax, Nr, a ,c)
    zeta = linspace(0,1,Nr);
    zeta = zeta';
    dzeta = zeta(2) - zeta(1);
    b = log(Rmax) - a - c;
    r = exp(a.*zeta.^3 + b.*zeta.^2 + c.*zeta);
    paran = 3*a.*zeta.^2 + 2*b.*zeta + c;
    d1 = 1./(paran.*r);
    d2 = - (6*a.*zeta + 2*b + paran.^2).*d1.^2./paran;
end

function [Bn_se, Bn_mo] = surface_velocity_expansion(feeding_area, mu2, Num_modes)
    
    M = 401;
    coo = linspace(-1,1,M);
    mu1 = 1-2*feeding_area; 
    right = mu1; 
    left = mu2; 
    
    [Bn0, Vn, u_design] = ExpandOnLegendre(Num_modes, 1, left, right);
    % sessile
    [~, energy_se] = sessileP(Bn0, Vn, coo, Num_modes);
    [Bn_se] = ExpandOnLegendre(Num_modes, 1/sqrt(energy_se), left, right);    
    % motile
    [~, energy_mo] = motileP(Bn0, Vn, coo, Num_modes);
    [Bn_mo] = ExpandOnLegendre(Num_modes, 1/sqrt(energy_mo), left, right);

end

function [us, NormP] = sessileP(Bn, Vn, coo, Num_modes)

    us = Bn(1).*Vn(1,coo);
    NormP = Bn(1)^2;% sessile
    for n = 2:Num_modes
        us = us + Bn(n).*Vn(n,coo);
        NormP = NormP + 2*Bn(n)^2/n/(n+1);
    end
    us(1) = 0;us(end) = 0;
end

function [us, NormP] = motileP(Bn, Vn, coo, Num_modes)
    us = Bn(1).*Vn(1,coo);
    NormP = 2*Bn(1)^2/3;% swim
    for n = 2:Num_modes
        us = us + Bn(n).*Vn(n,coo);
        NormP = NormP + 2*Bn(n)^2/n/(n+1);
    end
    us(1) = 0;us(end) = 0;
end


function [coefficients, Vn, f_design] = ExpandOnLegendre(numTerms, P_scale, low, high)
    % design a velocity shape
    x_scale = pi/(high-low); 
    f_design = @(x) P_scale.*sin(x_scale*(high-x));
    
    % Define the first derivative of Legendre polynomials
    legendreDerivative = @(n, x) n .* (x.*legendreP(n, x) - legendreP(n-1, x)) ./ (x.^2 - 1);
    Vn = @(n,x) 2/(n*(n+1)).*sqrt(1-x.^2).*legendreDerivative(n,x);
        
    % Initialize the coefficients vector
    coefficients = zeros(1, numTerms);
    
    % Calculate the coefficients
    for i = 1:numTerms

        % Integrate the product of sigmoid and the derivative of Legendre polynomial
        integ = integral(@(x) f_design(x) .* Vn(i, x), low, high);
        
        % Divide the integral by the norm of the Legendre polynomial
        coefficients(i) = integ * i*(i+1)*(2*i+1)/8;
    end
end

function [cc, Sh] = Advection_Diffusion_varycap_stretchmesh(r, d1, d2, dze, theta, ur, ut, Pe, indcap)
 
    Nr = length(r);
    Ntheta = length(theta);
    dr = r(2)-r(1);
    dtheta = theta(2)-theta(1);
    
    peur = ur.*Pe;
    peut = ut.*Pe;
    
    %- Initializing computation domain
    NumNodes = Nr * Ntheta;
    A = spalloc(NumNodes, NumNodes, 9*NumNodes);
    RHS = zeros(NumNodes,1);
    Node = zeros(Nr, Ntheta); 
    Node(1:NumNodes) = [1:NumNodes]; 
        

    for j = 1:Ntheta

        % Coeffecients

        ee(:,j) = -0.5*peur(:,j).*d1./dze - d1.^2./dze^2 + d1./r./dze + 0.5*d2./dze;
        gg(:,j) =  0.5*peur(:,j).*d1./dze - d1.^2./dze^2 - d1./r./dze - 0.5*d2./dze;
        ff(:,j) =  2.*d1.^2/dze^2 + 2/dtheta^2./r.^2 ;
        hh(:,j) = -0.5*peut(:,j)/dtheta./r - 1./r.^2/dtheta^2 + 0.5*cos(theta(j))/sin(theta(j))./r.^2/dtheta;
        pp(:,j) =  0.5*peut(:,j)/dtheta./r - 1./r.^2/dtheta^2 - 0.5*cos(theta(j))/sin(theta(j))./r.^2/dtheta;  
    end

    
    % Matrix for central points
    for i = 2:Nr-1
       for j = 2:Ntheta-1
          ANode_i = Node(i,j);
    
          A(ANode_i, Node(i-1, j  ) ) = ee(i,j);      
          A(ANode_i, Node(i  , j  ) ) = ff(i,j);
          A(ANode_i, Node(i+1, j  ) ) = gg(i,j);
          A(ANode_i, Node(i  , j-1) ) = hh(i,j);
          A(ANode_i, Node(i  , j+1) ) = pp(i,j);
       end
    end
        
    % boundary conditions-------------------------------------------------- 
    % boundary points
    ANode_i = Node(1,1);
    A(ANode_i, Node(1,1)) = 1;
    RHS(ANode_i) = 1;
    
    ANode_i = Node(Nr,1);
    A(ANode_i, Node(Nr,1)) = 1;
    RHS(ANode_i) = 0;
    
    ANode_i = Node(1,Ntheta); 
    A(ANode_i, Node(1,Ntheta)) = 1;
    ANode_i = Node(1,Ntheta-1); 
    A(ANode_i, Node(1,Ntheta-1)) = -1;
    RHS(ANode_i) = 0;
    
    ANode_i = Node(Nr,Ntheta);
    A(ANode_i, Node(Nr,Ntheta)) = 1;
    RHS(ANode_i) = 0;
    
    % Boundary (sphere surface and infinity distance)
    i = 1;
    for j = 2:indcap
        ANode_i = Node(i,j);
        A(ANode_i, Node(i,j)) = 1;
        RHS(ANode_i) = 1;
    end
    for j = indcap+1:Ntheta-1
        ANode_i = Node(i,j);
        A(ANode_i, Node(i,j)) = -1;
        A(ANode_i, Node(i+1,j)) = 1;
        RHS(ANode_i) = 0;
    end
    
    i = Nr;
    for j = 2:Ntheta
       ANode_i = Node(i,j);
       A(ANode_i, Node(i,j)) = 1;
       RHS(ANode_i) = 0;
    end
    
    % boundary (theta = 0 and theta = pi)
    j = 1;
    for i = 2:Nr-1
        ANode_i = Node(i,j); 
        A(ANode_i, Node(i,j)) = -1;
        A(ANode_i, Node(i,j+1)) = 1;
        RHS(ANode_i) = 0;
    end
    
    j = Ntheta;
    for i = 2:Nr-1
        ANode_i = Node(i,j);
        A(ANode_i, Node(i,j)) = 1;
        A(ANode_i, Node(i,j-1)) = -1;
        RHS(ANode_i) = 0;   
    end
    
    % concentration
    c = A\RHS;
    cc = reshape(c, Nr, Ntheta);

    surf_grad = sum(0.5*(-3.*cc(1,:)+4.*cc(2,:)-cc(3,:)).*sin(theta)./dr);
%     surf_grad = sum((-cc(1,:)+cc(2,:)).*sin(theta).*d1(1)./dze);
    Sh = -0.5*surf_grad*dtheta;

end
