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

#ifndef SPRINGBONDINFO_H
#define SPRINGBONDINFO_H

#include "BondInfo.h"

class SpringBondInfo:public BondInfo
{
	private:
		double k;
		double eqDist;
	public:
		SpringBondInfo();
		SpringBondInfo(double newK, double newEqDist);
		double calcForce(double dist);
        double getEqDist();
};

#endif