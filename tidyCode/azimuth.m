%% 计算方位角
function azi = azimuth( x1,y1,x2,y2 )

    if ((x2 >= x1) && (y2 == y1)) 
        azi=0;
        return
    end
    if ((x2 == x1) && (y2 >= y1)) 
        azi=0.5*pi;
        return
    end
    if ((x2 <= x1) && (y2 == y1)) 
        azi=pi;
        return
    end
    if ((x2 == x1) && (y2 <= y1)) 
        azi=1.5*pi;
        return
    end
    azi = (y1 - y2) / (x1 - x2);
    if ((x2 > x1) && (y2 > y1)) 
        azi = atan(azi);
    end
    if ((x2 <= x1) && (y2 > y1)) 
        azi = atan(azi) + pi;
    end
    if ((x2 <= x1) && (y2 <= y1)) 
        azi = atan(azi) + pi;
    end
    if ((x2 > x1) && (y2 <= y1)) 
        azi = atan(azi) + pi * 2;
    end

end