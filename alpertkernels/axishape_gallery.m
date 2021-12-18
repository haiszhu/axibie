function Z = axishape_gallery(shapeIn) 
%
% Hai 05/09/20

  ax = 1; az = 1;
  switch shapeIn
   case 'ellipseZ'
    rho = 1;
    ax = 1/2; 
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);
   
   case 'ellipseZ2'
    rho = 1;
    ax = 1/4; 
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);

   case 'ellipseX'
    rho = 1;
    az = .5;
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);
    
   case 'oblateZ'
    rho = 1;
    az = 1/3; 
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);

   case 'dumbbell'
    rho = @(t) 1 + real(Ynm(2,0,t(:),0*t(:)));
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));
    
%    case 'tilt_dumbbell'	
%     rho = @(t) 1 + real(Ynm(2,1,t(:),0*t(:))) + .1*real(Ynm(3,2,t(:),0*t(:)));
%     Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));
    
   case 'butternut'
    %rho = .1+abs(Ynm(1,0,u,v)).^2;
    cf = .8;
    s = @(t) 1-cf*cos(pi*(cos(t)-.2)/2); 
    Z = @(t) sqrt((1-cos(t).^2).*s(t).^2) - 1i*cos(t);
    
   case 'tri_ap'  % ~ triangle antiprism
    rho = @(t) 1 + .2*exp(-3*real(Ynm(3,2,t(:),0*t(:))));
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));
    
   case 'square_ap' % ~ square antiprism
    rho = @(t) exp(.5*sin(t).^4.*cos(t));
    Z = @(t) ax*rho(t).*sin(t) - 1i*az*rho(t).*cos(t);
    
   case 'oblate85'
    rho = 1;
    ax = .47;
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);
    
   case 'oblate75'
    rho = 1;
    ax = .366;
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);

   case 'oblate65'
    rho = 1;
    ax = .29;
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);
    
   case 'Y43'
    n = 4; m = 3; r=0.5;                      
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:)));      
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:)); 
    
   case 'Y87'
    n = 8; m = 7; r=0.5; 
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:)));         
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:)); 
    
    case 'Y1615'
    n = 16; m = 15; r=0.5;   
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:)));      
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:)); 
    
   case 'RY43'  
    n = 4; m = 3; 
    r = rand();                      
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:)));      
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));  
    
   case 'RY87'
    n = 8; m = 7;
    r = rand();
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:))); 
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));  
    
    case 'RY1615'
    n = 16; m = 15;
    r = rand();  
    rho = @(t) 1 + r*real(Ynm(n,m,t(:),0*t(:))); 
    Z = @(t) ax*rho(t).*sin(t(:)) - 1i*az*rho(t).*cos(t(:));  
   
    case 'stomatocyte'
    lam = 0.95; 
    Z = @(t) (-(1.5-cos(t)).*exp(-lam*1i*pi*sin(t)))/2.2;    
    
   otherwise
    rho = 1;
    Z = @(t) ax*rho.*sin(t) - 1i*az*rho.*cos(t);
    
  end
