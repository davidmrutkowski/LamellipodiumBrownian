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

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <omp.h>
#include <algorithm>
#include <set>
#include "SimulationBox.h"
#include "Grid.h"
#include "ParticleInfo.h"
#include "MiscStructs.h"
#include "Filament.h"
#include "SpringBondInfo.h"

using namespace std;

class Simulation
{
	private:
		const double pi = 3.14159265359;
        const double boltzmannConstant = 1.38e-23*1e12*1e6;
        
        double temperature;
        double dt;
        double zeta;
        double zetaInNascentFA;
        int numSteps;
        int numThreads;
        int snapshotStep;
        int updateGridStep;
        double filamentK;
        double persistenceLength;
        int numMonomersPerBead;
        double L_zero;
        double eta;
        double etaInNascentAdhesion;
        double criticalBendingForceMag;
        double criticalBondForceMag;
        double thermal_force;
        double thermal_force_inNascentFA;
        unsigned seed;
        double totalSimulationTime;
        double snapshotTime;
        double currSimTime;
        
        double delta0;
        double maxYValInitial;
        double initialDeltaOverDelta0;
        double finalDeltaOverDelta0;
        double deltaOverDelta0StepDown;
        int numStepsBetweenDeltaSteps;
        double omega;
        double filamentBondK;
        double branchAngleK;
        double wallAngleK;
        double pullForceMag;
        
        double minCrosslinkBondDist;
        double maxCrosslinkBondDist;
        
        Coordinate FABeadPos;
        double FABeadRadius;
        double FABeadLength;
        
        double rRepulsiveInteraction;
        
        bool thermalForcesOn;
        
        bool severCrosslinksWithExtension;
        bool crosslinkNeighboringUponAdding;
        bool filamentBondsInBondList;
        bool excludedVolOn;
        
        bool uniformPullingForceOn;
        bool backPullingForceOn;
        
        std::vector <int> removalList;
        
        std::vector <std::vector <int>> excludedVolumeList;
        
        std::string outputFileName;
        
        std::string positionFileName;
        std::string bondFileName;
        std::string angleFileName;
        
        SimulationBox simBox;
        Grid particleGrid;
        
        ParticleInfo pinfo;
        
        std::vector <Filament> filaments;
        std::vector <SpringBondInfo> bondTypes;
        
        TaggedVector f_TaggedVector;
        
        std::uniform_real_distribution<> dis;
        std::normal_distribution<> gaussianDis;
    
        std::vector <std::mt19937> generators;
        
        std::vector<std::vector <Coordinate>> tempForces;
        
        std::vector<std::vector <Coordinate>> tempStresses;
        
        std::deque <int> bondTagsOldestToYoungest;
        
        double nascentFABindProbability;
        double nascentFAUnbindProbability;

	public:
		Simulation();
		
		void readParameterFile(std::string file_name);
        void run();
        void writeXYZFile(std::ofstream& outputXYZ);
        void writeBondFile(std::ofstream& outputBND);
        double calcForces();
        void calcFilamentForces(int f);
        void calcBondForce(Bond b);
        Coordinate getBondForce(Bond b);
        Coordinate getBondForce(Coordinate ij, int bondType);
        void calcFilamentAngleForce(Angle ang);
        void calcAngleForce(Angle ang);
        Coordinate calcExcludedVolumeForces(Coordinate& ci, Coordinate& cj, Coordinate& ci2, Coordinate& cj2);
        void putAllParticlesInGrid();
        void moveParticles();
        void readAngleFile();
        void readBondFile();
        void readXYZ();
        void initDerivedValues();
        
        void removeParticlesOutsideBox();
        void removeParticleByIndex(int pIndex);
        void removeFilament(int filamentTag);
        void removeFilamentsNearPos(Coordinate pos, double dist);
        void applyNascentFAForce(Coordinate pos, double dist);
        void removeAllCrosslinks();
        void removeTemporaryCrosslinks();
        void removeBond(int bondTag);
        void removeBondByIndex(int bondIndex);
        int removeFlaggedBonds();
        int addTemporaryCrosslinks(int nCrosslinks);
        
        void addRandomFilament(bool crosslinkNeighboring);
        void addFilament(int numBeads, Coordinate startPos, Coordinate direction, bool crosslinkNeighboring);
        int addParticlePutInGrid(Particle p, int numStressMeasurements);
        double calcLeadingEdgeForce(int f);
        int addCrosslinksInOrder(int nCrosslinks);
        
        void cleanExistingFilament(int filamentIndex, int filamentTag);
        void cleanNewFilament(Filament tempFil);
        
        void applyNascentFAForceToTags(std::vector <int> tagList);
        
        double calcMinDistanceCylinders(Coordinate ri, Coordinate rj, Coordinate ri2, Coordinate rj2);
};


#endif