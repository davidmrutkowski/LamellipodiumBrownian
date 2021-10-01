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

#ifndef SIMULATIONBOX_H
#define SIMULATIONBOX_H

#include <iostream>
#include <stdexcept>
#include "MiscStructs.h"

using namespace std;

class SimulationBox
{
	private:
		Coordinate boxLengths;
        bool periodicX, periodicY, periodicZ;
        
	public:
        SimulationBox();
		SimulationBox(Coordinate b, bool px, bool py, bool pz);
		
        double periodicWrap(double pos, double boxL);
        Coordinate periodicWrap(Coordinate a);
        struct Coordinate calcDisplacement(Coordinate c, Coordinate d);
        double calcDistance(Coordinate c, Coordinate d);
        double getBoxLength(int index);
        void setBoxLength(int index, double newVal);
};

#endif