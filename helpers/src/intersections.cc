#include "intersections.hh"

#include <iostream>
#include <stdexcept>

namespace aidaTT
{
    circle::circle(double c0, double c1, double r) : _c0(c0), _c1(c1), _r(r)
    {
        if(r < 0.)
            throw std::invalid_argument("[circle::circle] Error: radius must be >0  when creating a circle.");
    }



    straightLine::straightLine(double n0, double n1, double d) : _n0(n0), _n1(n1), _d0(d)
    {
        if(d < 0.)
            throw std::invalid_argument("[straightLine::straightLine] Error: distance must be >= 0 when creating a straight line.");

        if(fabs(n0) < 1e-12 && fabs(n1) < 1e-12)
            {
                throw std::invalid_argument("[straightLine::straightLine] Error: normal vector must be != zero.");
            }
    }



    intersections intersectCircleCircle(const circle& circ1, const circle& circ2)
    {
        double x1 = circ1.center().first;
        double y1 = circ1.center().second;
        double r1 = circ1.radius();

        double x2 = circ2.center().first;
        double y2 = circ2.center().second;
        double r2 = circ2.radius();

        /// reduce the two circles to a straight line
        double nx = 2 * (x2 - x1);
        double ny = 2 * (y2 - y1);
        double dist = r1 * r1 - r2 * r2 - x1 * x1 + x2 * x2 - y1 * y1 + y2 * y2;

        straightLine reducedCircle(nx, ny, dist);
        return intersectCircleStraightLine(circ1, reducedCircle);
    }



    intersections intersectCircleStraightLine(const circle& c,  const straightLine& line)
    {
        const double x0 = c.center().first;
        const double y0 = c.center().second;
        /// adjust straight line to use the same coordinates
        straightLine sL(line);
        sL.move(x0, y0);

        // check if it does not intersect or just touch
        const double discriminant = c.r2() * sL.normalSquare() - sL.d2();

        if(discriminant  < 0. &&  fabs(discriminant) > 1e-9)
            return intersections();

        const double nx = sL.normal().first;
        const double ny = sL.normal().second;
        const double d = sL.distance();
        const double DISC = sqrt(discriminant);

        const double x1 = x0  + (nx * d + ny * DISC) / sL.normalSquare();
        const double y1 = y0  + (ny * d - nx * DISC) / sL.normalSquare();

        intersections retType;
        retType.add(x1, y1);

        /// return only a single solution, if both are very close
        if(fabs(discriminant) <= 1e-9)
            return retType;
        const double x2 = x0  + (nx * d - ny * DISC) / sL.normalSquare();
        const double y2 = y0  + (ny * d + nx * DISC) / sL.normalSquare();
        retType.add(x2, y2);
        return retType;
    }



    intersections intersectStraightLineStraightLine(const straightLine& line1, const straightLine& line2)
    {
        const double nx1 = line1.normal().first;
        const double ny1 = line1.normal().second;
        const double d1 = line1.distance();

        const double nx2 = line2.normal().first;
        const double ny2 = line2.normal().second;
        const double d2 = line2.distance();

        intersections retType;

        /// test linear independency
        const double discriminant = nx1 * ny2 - nx2 * ny1;
        if(discriminant < 1e-9)
            return retType;

        const double x = (d1 * ny2 - ny1 * d2) / discriminant;
        const double y = (d2 * nx1 - d1 * nx2) / discriminant;
        retType.add(x, y);

        return retType;
    }
}
