#include "intersections.hh"

#include <iostream>
namespace aidaTT
{
    intersections intersectCircleCircle(const circle& circ1, const circle& circ2)
    {
        return intersections();

    }



    intersections intersectCircleStraightLine(const circle& c,  const straightLine& line)
    {
        const double x0 = c.center().first;
        const double y0 = c.center().second;
        /// adjust straight line to use the same coordinates
        straightLine sL(line);
        sL.move(x0, y0);

        // check if it does not intersect or just touch
        const double discriminant = c.r2() - sL.d2();

        if(discriminant  < 0. &&  fabs(discriminant) > 1e-6)
            return intersections();

        const double nx = sL.normal().first;
        const double ny = sL.normal().second;
        const double d = sL.distance();
        const double DISC = sqrt(discriminant);

        const double x1 = x0  + (nx * d + ny * DISC);
        const double y1 = y0  + (ny * d - nx * DISC);

        intersections retType;
        retType.add(x1, y1);

        /// return only a single solution, if both are very close
        if(fabs(discriminant) < 1e-9)
            return retType;
        const double x2 = x0  + (nx * d - ny * DISC);
        const double y2 = y0  + (ny * d + nx * DISC);
        retType.add(x2, y2);
        return retType;
    }



    intersections intersectStraightLineStraightLine(const straightLine& line1, const straightLine& line2)
    {
        return intersections();
    }
}
