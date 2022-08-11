function Elliptic_Confidence_Region(a,b,cx,cy,angle,color,mode)

for m = 2:length(angle)
    if angle(m)-angle(m-1) > pi
        angle(m) = angle(m) - 2*pi;
    end
    if angle(m)-angle(m-1) < -pi
        angle(m) = angle(m) + 2*pi;
    end
end

CX = linspace(cx(1),cx(end),1000)';
A = interp1(cx,a,CX);
B = interp1(cx,b,CX);
CY = interp1(cx,cy,CX);
ANGLE = interp1(cx,angle,CX);

RRR = 100;
Points = zeros(1000*RRR,2);

for m = 1:1000
    r = linspace(0,2*pi+0.1,RRR);
    
    p = [(A(m)*cos(r))',(B(m)*sin(r))'];
    
    beta = [cos(ANGLE(m)),sin(ANGLE(m));-sin(ANGLE(m)),cos(ANGLE(m))];
    
    Points((m-1)*RRR+1:RRR*m,:) = p*beta;
    Points((m-1)*RRR+1:RRR*m,1) = Points((m-1)*RRR+1:RRR*m,1) + CX(m);
    Points((m-1)*RRR+1:RRR*m,2) = Points((m-1)*RRR+1:RRR*m,2) + CY(m);
end

K = boundary(Points(:,1),Points(:,2),1);
if strcmp(mode,'sed_rate') == 1
    patch(Points(K,1),exp(Points(K,2)),1,'FaceColor',color,'EdgeColor','none');
elseif strcmp(mode,'log_sed_rate') == 1 || strcmp(mode,'lag') == 1
    patch(Points(K,1),Points(K,2),1,'FaceColor',color,'EdgeColor','none');
end
alpha(.2);


end

