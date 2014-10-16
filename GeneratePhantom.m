function [ ImageOut ] = GeneratePhantom(  )
%Generates tubular phantom with artificial break points.

ImageOut = zeros(20,20,50);

for z = 5:20
    for x = 5:15
        for y = 5:15
            if (x-10)^2 + (y-10)^2 < 15
                ImageOut(z,x,y) = 1;
            end
        end
    end
end

for z = 25:45
    for x = 5:15
        for y = 5:15
            if (x-10)^2 + (y-10)^2 < 15
                ImageOut(z,x,y) = 1;
            end
        end
    end
end




end

