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

#include <iostream>
#include "Particle.h"

using namespace std;

Particle::Particle()
{
	int dimensions = 3;
	
	for(int i = 0; i < dimensions; i++)
	{
		r[i] = 0.0;
		force[i] = 0.0;
	}
	
	mass = 0.0;
    bin = -1;
    type = -1;
}

Particle::Particle(const Particle &p2)
{
	for(int i = 0; i < 3; i++)
	{
		r[i] = p2.r[i];
		force[i] = p2.force[i];
	}
	
	mass = p2.mass;
    bin = p2.bin;
    type = p2.type;
}

Particle::Particle(double rx, double ry, double rz)
{
	r[0] = rx;
	r[1] = ry;
	r[2] = rz;
	
	int dimensions = 3;
	
	for(int i = 0; i < dimensions; i++)
	{
		force[i] = 0.0;
	}
	
	mass = 0.0;
    bin = -1;
    type = -1;
}

void Particle::setPosition(int index, double value)
{
	if(index < 3)
		r[index] = value;
}

double Particle::getPosition(int index)
{
	if(index < 3)
		return r[index];
	else
		return -1.0;
}

void Particle::setForce(int index, double value)
{
	if(index < 3)
		force[index] = value;
}

void Particle::addForce(int index, double value)
{
	if(index < 3)
		force[index] += value;
}

void Particle::zeroAllForces()
{
	for(int i = 0; i < 3; i++)
	{
		force[i] = 0.0;
	}
}

double Particle::getForce(int index)
{
	if(index < 3)
		return force[index];
	else
		// maybe should throw exception here
		return 0.0;
}

double Particle::getMass()
{
	return mass;
}

void Particle::setBin(int newBin)
{
	bin = newBin;
}

int Particle::getBin()
{
	return bin;
}

void Particle::setType(int newType)
{
    type = newType;
}
int Particle::getType()
{
    return type;
}