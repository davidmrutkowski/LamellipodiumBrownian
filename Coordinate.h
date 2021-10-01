/**
 * Copyright (C) 2021 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
 */
 
#ifndef COORDINATE_H
#define COORDINATE_H

#include <math.h>

struct Coordinate
{
	double x, y, z;
    
    Coordinate()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
    
    Coordinate (double a, double b, double c) : x(a), y(b), z(c)
    {
    }
    
	struct Coordinate getUnitCoord()
	{
		double mag = sqrt(x*x + y*y + z*z);
		
		double invMag = 1.0/mag;
		
		struct Coordinate uCoord = {x*invMag, y*invMag, z*invMag};
		
		return uCoord;
	}
	
	double getMagnitude()
	{
		return sqrt(x*x + y*y + z*z);
	}
	
	Coordinate crossProduct(Coordinate b)
	{
		Coordinate c;
		
		c.x = this->y * b.z - this->z * b.y;
		c.y = this->z * b.x - this->x * b.z;
		c.z = this->x * b.y - this->y * b.x;
		
		return c;
	}
	
	
	Coordinate operator+(const Coordinate& b)
	{
		Coordinate c;
		c.x = this->x + b.x;
		c.y = this->y + b.y;
		c.z = this->z + b.z;
		return c;
	}
	
	Coordinate operator-(const Coordinate& b)
	{
		Coordinate c;
		c.x = this->x - b.x;
		c.y = this->y - b.y;
		c.z = this->z - b.z;
		return c;
	}
	
	Coordinate operator-() const
	{
		Coordinate c;
		
		c.x = -this->x;
		c.y = -this->y;
		c.z = -this->z;
		return c;
	}
	
	double operator*(Coordinate b) const
	{
		return this->x * b.x + this->y * b.y + this->z * b.z;
	}
	
	Coordinate operator*(double a) const
	{
		Coordinate c;
		c.x = this->x * a;
		c.y = this->y * a;
		c.z = this->z * a;
		return c;
	}
	
	Coordinate operator/(double a) const
	{
		Coordinate c;
		c.x = this->x / a;
		c.y = this->y / a;
		c.z = this->z / a;
		return c;
	}
};

inline Coordinate operator*(const double a, const Coordinate& b)
{
	Coordinate c = b * a;
	return c;
}

#endif