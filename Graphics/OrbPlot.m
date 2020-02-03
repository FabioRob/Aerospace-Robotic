function [] = OrbPlot(r0,v0,mu,col,linewidth,TA_final)

% -----------------------------------------------------------------------%
%
% OrbPlot plots the orbit from state vector
%
% Arguments:
%
% r0        - position vector at t0
% v0        - velocity vector at t0
% mu        - gravitational paramater of the central body(Km^3/s^2)
% col       - line color spec
% linewidth - line width spec
% TA_final -  indicates orbit's portion to be plotted
%             if TA_final is not provided, it plots the entire orbit
%
% -----------------------------------------------------------------------%

deg=pi/180;

R0=norm(r0);
vr0=dot(r0,v0)/R0;
coe = coe_from_sv(r0,v0,mu);
H=coe(1);
e = coe(2);
TA_in=(zero_to_360( coe(6)/deg ) )*deg ; % true anomaly at departure in [0-2pi](rad)
p=H^2/mu;
len=0;
if nargin==5
    if e<1
        dT=(0:0.000001:2*pi);         
    else % hyperbolic orbit(only from/to periapsis)
        TA_inf = acos(-1/e);
        dT=(0:0.000001:(TA_inf-0.5));
    end
    
elseif e<1
    TA_final=(zero_to_360( TA_final/deg ) )*deg ; % true anomaly at arrival in [0-2pi](rad)
    if TA_final<TA_in
     dT=(0 :0.000001:abs(2*pi-TA_in) ) ;  
     len=length(dT);
     dT = [dT , ((TA_final-TA_in) : -0.000001: -TA_in) ];
    else
     dT=(0 :0.000001:abs(TA_final-TA_in) ) ;     
    end
else % hyperbolic orbit(only from/to periapsis)   
     dT=(0: 0.000001*sign(TA_final) :TA_final);   
end

rdeltaT=p./(1+cos(dT).*(H^2/(mu*R0)-1)-sin(dT).*(H*vr0/mu));  % norm of position vector in function of variation in the true anomaly from r0 (eq. 2.152 Curtis)
f=1-(rdeltaT./p).*(1-cos(dT)); % Lagrange coefficient f
g=((rdeltaT.*R0)./H).*sin(dT); % Lagrange coefficient g
 % position vectors for all anomalies in the given range
rTX=f.*r0(1,1)+g.*v0(1,1);
rTY=f.*r0(end-1)+g.*v0(end-1);
rTZ=f.*r0(end)+g.*v0(end);
hold all
if len==0
    plot3(rTX,rTY,rTZ,'color',col,'linewidth',linewidth)
else 
    plot3(rTX(1:len),rTY(1:len),rTZ(1:len),'color',col,'linewidth',linewidth)
    plot3(rTX(len+1:end),rTY(len+1:end),rTZ(len+1:end),'color',col,'linewidth',linewidth)
end

xlabel('X')
ylabel('Y')
zlabel('Z')
end

