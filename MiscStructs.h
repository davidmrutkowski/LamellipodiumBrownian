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

#ifndef MISCSTRUCTS_H
#define MISCSTRUCTS_H

#include "Coordinate.h"

struct Bond
{
	int i, j, bondTypeIndex;
    double creationTime;
    double destructionTime;
    
    Bond(int a, int b, int c) : i(a), j(b), bondTypeIndex(c)
    {
        creationTime = 0.0;
    }
    
    Bond(int a, int b, int c, double d) : i(a), j(b), bondTypeIndex(c), creationTime(d)
    {
    }
    
    Bond(int a, int b, int c, double d, double e) : i(a), j(b), bondTypeIndex(c), creationTime(d), destructionTime(e)
    {
    }
    
    bool operator<(const Bond& a) const
    {
        return (this->i < a.i);
    }
};

struct HalfBond
{
    int j, bondTag;
    
    HalfBond(int a, int b) : j{a}, bondTag{b}
    {
    }
};

struct Angle
{
	int i, j, k, angleTypeIndex;
    
    Angle(int a, int b, int c, int d) : i{a}, j{b}, k{c}, angleTypeIndex{d}
    {
    }
};

struct TemporaryBond
{
	int i, j;
	struct Coordinate force;
};

struct TemporaryBondAndDistance
{
    int i, j;
    struct Coordinate force;
    struct Coordinate distance;
};

struct TemporaryBondMag
{
    int i, j;
    double forceMag;
};

struct TemporaryForce
{
	int i;
	struct Coordinate force;
};

#endif