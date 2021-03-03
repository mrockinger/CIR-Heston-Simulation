function q = inquadrant(z)
% returns the quadrant in which z is located
theta=angle(z);
if theta>0 && theta<=pi/2
    q=1;
elseif theta>pi/2 && theta<=pi
    q=2;
elseif theta<=0 && theta>=-pi/2
    q=4;
elseif theta<-pi/2 && theta>=-pi
    q=3;
end
    