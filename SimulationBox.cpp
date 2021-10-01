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

#include "SimulationBox.h"

SimulationBox::SimulationBox()
{
    boxLengths = Coordinate {1.0, 1.0, 1.0};
    periodicX = true;
    periodicY = true;
    periodicZ = true;
}

SimulationBox::SimulationBox(Coordinate b, bool px, bool py, bool pz)
{
    boxLengths = b;
    
    periodicX = px;
    periodicY = py;
    periodicZ = pz;
}

double SimulationBox::getBoxLength(int index)
{
    if(index == 0)
        return boxLengths.x;
    else if(index == 1)
        return boxLengths.y;
    else if(index == 2)
        return boxLengths.z;
    else
        throw std::runtime_error("SimulationBox does not have an index of " + std::to_string(index));
    
    return 0.0;
}

void SimulationBox::setBoxLength(int index, double newVal)
{
    if(index == 0)
        boxLengths.x = newVal;
    else if(index == 1)
        boxLengths.y = newVal;
    else if(index == 2)
        boxLengths.z = newVal;
}

double SimulationBox::periodicWrap(double pos, double boxl)
{
    return pos - boxl * round(pos/boxl);
}

Coordinate SimulationBox::periodicWrap(Coordinate a)
{
    if(periodicX)
        a.x = periodicWrap(a.x, boxLengths.x);
    if(periodicY)
        a.y = periodicWrap(a.y, boxLengths.y);
    if(periodicZ)
        a.z = periodicWrap(a.z, boxLengths.z);
    
    return a;
}

struct Coordinate SimulationBox::calcDisplacement(Coordinate a, Coordinate b)
{
    Coordinate c = a - b;
    
    c = periodicWrap(c);
    
    return c;
}

double SimulationBox::calcDistance(Coordinate a, Coordinate b)
{
    Coordinate c = calcDisplacement(a, b);
    
    return c.getMagnitude();
}