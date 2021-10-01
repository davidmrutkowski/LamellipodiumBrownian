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
 
#ifndef FILAMENT_H
#define FILAMENT_H

#include <vector>
#include <queue>
#include <iostream>
#include "Coordinate.h"
#include "ParticleInfo.h"

class Filament
{
	private:
		std::deque <int> p_tags;
		bool isDaughter;
		int motherTag;
        bool cappedBarbedEnd;
        bool cappedPointedEnd;
	public:
		Filament();
		
		int getNumParticles();
		int getTagAtIndex(int pos);
		
		int addParticleFront(int newTag);
		int addParticleBack(int newTag);
		void removeParticleFront();
		void removeParticleBack();
		
		struct Coordinate calcBendingPotential(double k, ParticleInfo &pinfo, struct TemporaryBondMag *maxBond, bool bendingInStress);

		struct Coordinate getFrontDirection(ParticleInfo pinfo);
		
		static struct Coordinate matrixProduct(double basis[3][3], struct Coordinate a);
		
        void setMotherTag(int newMotherTag);
		int getMotherTag();
        
		void setIsDaughter(bool state);
		bool getIsDaughter();
		
        Filament breakFilamentAtPos(int pos);
        Filament removeBeadAtPos(int pos);
        Filament removeBeadByTag(int tag);
        
        int getPointedTag();
        int getBarbedTag();
        
        void setCappedBarbedEnd(bool newState);
        bool getCappedBarbedEnd();
        void setCappedPointedEnd(bool newState);
        bool getCappedPointedEnd();
        
        int getIndexOfTag(int tempTag);
        
        std::deque <int> getAllTags();
};

#endif