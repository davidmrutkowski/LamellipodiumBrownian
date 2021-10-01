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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <map>
#include <queue>
#include <array>
#include <math.h>
#include "Coordinate.h"

class Particle
{
	private:
		std::array <double, 3> r;
		std::array <double, 3> force;
		double mass;
		int bin;
        int type;
		
	public:
		Particle();
		Particle(double rx, double ry, double rz);
		Particle(const Particle &p2);
		void setPosition(int index, double value);
		double getPosition(int index);
		void setForce(int index, double value);
		double getForce(int index);
		void zeroAllForces();
		void addForce(int index, double value);
		double getMass();
		void setBin(int newBin);
		int getBin();
        void setType(int newType);
        int getType();
};

#endif // USER_H