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

#include "Simulation.h"

Simulation::Simulation()
{
    // default values
    temperature = 300.0;
    eta = 0.301;
    etaInNascentAdhesion = eta;
    numMonomersPerBead = 0;
    dt = 0.1;
    totalSimulationTime = 1.0;
    snapshotTime = 0.01;
    persistenceLength = 17.0;
    numThreads = 1;
    thermalForcesOn = true;
    branchAngleK = 0.1;
    filamentBondK = 100.0;
    //pullForceMag = 0.00135;
    pullForceMag = 0.00135*0.5*3;
    // diameter of rod
    rRepulsiveInteraction = 7e-3;
    severCrosslinksWithExtension = false;
    crosslinkNeighboringUponAdding = true;
    filamentBondsInBondList = true;
    excludedVolOn = false;
    backPullingForceOn = true;
    uniformPullingForceOn = false;
    
    FABeadPos = {0.0, 0.0, -1.0-0.125-0.5};
    FABeadRadius = 0.25*0.5;
    FABeadLength = 1.25 - FABeadRadius * 2;
    
    nascentFABindProbability = 0.5;
    nascentFAUnbindProbability = 0.5;
    
    
    currSimTime = 0.0;
    
    outputFileName = "";
    
    initDerivedValues();

    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	// this does not randomize the seed on windows!
    //std::random_device rd;
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    dis = std::uniform_real_distribution<> (0.0, 1.0);
    gaussianDis = std::normal_distribution<> (0.0, 1.0);
    
    //https://stackoverflow.com/questions/37305104/how-to-generate-random-numbers-in-a-thread-safe-way-with-openmp
	for (int i = 0; i < omp_get_max_threads(); i++)
	{
		generators.emplace_back(std::mt19937(gen));
	}
    
    //BondInfo binfo;
    //AngleInfo ainfo;
}

void Simulation::initDerivedValues()
{
    L_zero = 2.7e-3 * numMonomersPerBead;
    zeta = 4*pi*eta*L_zero / (0.84 + log(L_zero/2/0.0035));
    snapshotStep = (int)round(snapshotTime / dt);
    numSteps = (int)round(totalSimulationTime / dt);
    filamentK = temperature * boltzmannConstant * persistenceLength;
    thermal_force = sqrt(2*boltzmannConstant*temperature*zeta/dt);
    zetaInNascentFA = 4*pi*etaInNascentAdhesion*L_zero / (0.84 + log(L_zero/2/0.0035));
    thermal_force_inNascentFA = sqrt(2*boltzmannConstant*temperature*zetaInNascentFA/dt);
}

void Simulation::readParameterFile(std::string fileName)
{
    ifstream inputfile (fileName);
    double timeBetweenDeltaSteps = 0.0;
    
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        std::string line;
              
		while(std::getline(inputfile, line))
		{	
			std::istringstream iss(line);
			
			std::string category;
            std::string value;

            iss >> category >> value;

            // this should probably be a switch-case but cannot do on std::string
            if(category.compare("dt") == 0)
            {
                dt = std::stod(value);
            }
            else if(category.compare("totalSimulationTime") == 0)
            {
                totalSimulationTime = std::stod(value);
            }
            else if(category.compare("snapshotTime") == 0)
            {
                snapshotTime = std::stod(value);
            }
            else if(category.compare("persistenceLength") == 0)
            {
                persistenceLength = std::stod(value);
            }
            else if(category.compare("numMonomersPerBead") == 0)
            {
                numMonomersPerBead = std::stod(value);
            }
            else if(category.compare("temperature") == 0)
            {
                temperature = std::stod(value);
            }
            else if(category.compare("eta") == 0)
            {
                eta = std::stod(value);
            }
            else if(category.compare("positions") == 0)
            {
                positionFileName = value;
            }
            else if(category.compare("bonds") == 0)
            {
                bondFileName = value;
            }
            else if(category.compare("angles") == 0)
            {
                angleFileName = value;
            }
            else if(category.compare("numThreads") == 0)
            {
                numThreads = std::stod(value);
            }
            else if(category.compare("thermalForcesOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> thermalForcesOn;
            }
            else if(category.compare("filamentBondK") == 0)
            {
                filamentBondK = std::stod(value);
            }
            else if(category.compare("branchAngleK") == 0)
            {
                branchAngleK = std::stod(value);
            }
            else if(category.compare("severCrosslinksWithExtension") == 0)
            {
                istringstream(value) >> std::boolalpha >> severCrosslinksWithExtension;
            }
            else if(category.compare("outputName") == 0)
            {
                outputFileName = value;
            }
            else if(category.compare("excludedVolOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> excludedVolOn;
            }
            else if(category.compare("backPullingForceOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> backPullingForceOn;
            }
            else if(category.compare("uniformPullingForceOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> uniformPullingForceOn;
            }
            else if(category.compare("etaInNascentAdhesion") == 0)
            {
                etaInNascentAdhesion = std::stod(value);
            }
            else if(category.compare("nascentFABindProbability") == 0)
            {
                nascentFABindProbability = std::stoi(value);
            }
            else if(category.compare("nascentFAUnbindProbability") == 0)
            {
                nascentFAUnbindProbability = std::stoi(value);
            }
            else if(category.compare("reactionRateNascentFABind") == 0)
            {
                double reactionRate = std::stoi(value);
                nascentFABindProbability = 1.0 - exp(-reactionRate * dt);
            }
            else if(category.compare("reactionRateNascentFAUnbind") == 0)
            {
                double reactionRate = std::stoi(value);
                nascentFAUnbindProbability = 1.0 - exp(-reactionRate * dt);
            }
            else
            {
                throw std::runtime_error("Unknown value in input file associated with: " + category);
            }
            
            linecount++;
			
		}
        
		inputfile.close();
    }
    else
    {
        // throw an error
        throw std::runtime_error("main could not open: " + fileName);
    }
    
    initDerivedValues();
    
    Coordinate boxSize = {2.0, 0.3, 10.0};
    bool periodicX = true;
    bool periodicY = false;
    bool periodicZ = false;
    
    simBox = SimulationBox(boxSize, periodicX, periodicY, periodicZ);
    particleGrid = Grid(boxSize, periodicX, periodicY, periodicZ);
    particleGrid.setBinSize(L_zero*1.2);
    
    if(!positionFileName.empty())
    {
        readXYZ();
    }
    
    if(bondFileName.empty() && !positionFileName.empty())
    {
        
        // assumes that bond and angle file names are the same as positionFileName just with a different extension
        bondFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".bnd";
        
        cout << "Assuming bond file is named: " << bondFileName << endl;
    }
    readBondFile();
    
    if(angleFileName.empty() && !positionFileName.empty())
    {
        angleFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".ang";
        
        cout << "Assuming angle file is named: " << angleFileName << endl;
    }
    readAngleFile();
    
    /*if(filaments.size() > 0 && filaments[0].getNumParticles() > 1)
    {
        double lengthBetweenFirstAndSecond = (pinfo.getPos(filaments[0].getPointedTag()) - pinfo.getPos(filaments[0].getTagAtIndex(1))).getMagnitude();
        
        int estimatedNumMonomersPerBead = (int)round(lengthBetweenFirstAndSecond / 2.7e-3);
        
        if(numMonomersPerBead == 0)
        {
            numMonomersPerBead = estimatedNumMonomersPerBead;
            
            cout << "Assuming numMonomersPerBead as: " << numMonomersPerBead << endl;
        }
        else if(numMonomersPerBead != estimatedNumMonomersPerBead)
        {
            cout << "ERROR, MAY HAVE ISSUES: Estimated numMonomersPerBead is " << estimatedNumMonomersPerBead << " but inputed " << numMonomersPerBead << endl;
        }
    }*/
    
    // should determine this value by how fast individual particles are moving!!!
    double updateGridTime = 0.05;
    updateGridStep = (int)round(updateGridTime / dt);
    
    criticalBendingForceMag = 2.0;
    criticalBondForceMag = -2.0;
    
    minCrosslinkBondDist = 0.030;
    maxCrosslinkBondDist = 0.040;
    
    // filament bonds
    bondTypes.emplace_back(SpringBondInfo(filamentBondK, L_zero));
    
    // permanent crosslink bonds
    bondTypes.emplace_back(SpringBondInfo(100.0, 0.5*(minCrosslinkBondDist + maxCrosslinkBondDist)));
    
    // temporary crosslink bonds
    bondTypes.emplace_back(SpringBondInfo(100.0, 0.5*(minCrosslinkBondDist + maxCrosslinkBondDist)));
    
    thermal_force = sqrt(2*boltzmannConstant*temperature*zeta/dt);
}

void Simulation::readXYZ()
{
    ifstream inputfile (positionFileName);
    
    int numParticles = 0;
    
    if (inputfile.is_open())
	{
		int linecount = 0;
        
        int currFilamentType = -1;
        int newFilamentTag = -1;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numParticles;
            }
            else if(linecount == 1)
            {
                std::string tempString;
                iss >> tempString;
                
                int pos = tempString.find("=");
                
                std::string timeString = tempString.substr(pos+1,tempString.length());
                
                currSimTime = std::stod(timeString);
                //exit(0);
            }
            else if(linecount > 1)
            {
                int tag, type;
                double x,y,z;
                int viscosityType;
                
                iss >> type >> tag >> x >> y >> z >> viscosityType;
                
                if(type != currFilamentType)
				{
                    // not the same filament
                    // assumes beads in the same filament have the same type
					Filament newFil;
					newFil.setIsDaughter(false);
					filaments.push_back(newFil);
                    newFilamentTag = f_TaggedVector.add();
					
                    currFilamentType = type;
				}

				Coordinate newCoordinate = {x, y, z};
				
				Particle newParticle(x, y, z);
                // set type of particle to the filament tag
                // filament tag does not currently update so that they are conserved
                // but this is probably ok since only ever check that beads are on the same filament or not
                newParticle.setType(newFilamentTag);
                
                int numStressMeasurements = 2;
                
				int newTag = pinfo.addParticle(newParticle, numStressMeasurements);
                if(newTag != tag)
                {
                    // set tag of this particle to its correct value listed in the xyz file
                    pinfo.setTagAtPos(pinfo.getNumParticles()-1, tag);
                    newTag = tag;
                }
                
                pinfo.setViscosityTypeByIndex(pinfo.getNumParticles()-1, viscosityType);
                
				int priorTag = filaments[newFilamentTag].addParticleBack(newTag);
                    
                // add bonds here for filament as well? (needed to keep track of excluded volume interactions)
                if(!filamentBondsInBondList && priorTag > -1)
				{
					// adds a bond between priorTag and newTag of type 0 if they are on the same filament
					// as indicated from the type column in the xyz file
					Bond newBond = {priorTag, newTag, 0};
					
					pinfo.addBond(newBond);
				}
            }
            
            linecount++;
        }
        
        if(linecount-2 != numParticles)
        {
            // warning that the number of particles in the simulation is different from numParticles
            std::cout << "Error. " << numParticles << " particles desired in simulation but " << linecount-2 << " particles in " << this->positionFileName;
        }
        
        inputfile.close();
    }
    else
    {
        cout << "Could not open " + positionFileName << endl;
    }
}

void Simulation::readBondFile()
{
    ifstream inputfile (bondFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numBonds = 0;
        int countBonds = 0;
        
        int numParticles = pinfo.getNumParticles();
        
        std::string line;
        
        int maxCurrTag = pinfo.getCurrMaxTag();
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numBonds;
            }
            else if(linecount > 1)
            {
                int particleI, particleJ, bondType;
                double bond_creationTime, bond_destructionTime;
                
                iss >> bondType >> particleI >> particleJ >> bond_creationTime >> bond_destructionTime;
                
                if(particleI < 0 || particleI > maxCurrTag || particleJ < 0 || particleJ > maxCurrTag)
                {
                    // warning about particleI and/or particleJ being out of bounds
                    throw std::runtime_error(std::to_string(particleI) + " or " + std::to_string(particleJ) + " tags are larger than current max tag " + std::to_string(maxCurrTag));
                }
                
                //add bond to binfo
                //int bondTag = binfo.addBond(particleI, particleJ, bondType);
                
                // SHOULD SET BOND CREATION TIME HERE IF IT IS IN BND FILE!!!
                
                //add bond tag to both particleI and particleJ
                pinfo.addBond(Bond {particleI, particleJ, bondType, bond_creationTime, bond_destructionTime});
                countBonds++;
            }
            
            linecount++;
        }

        if(countBonds != numBonds)
        {
            // warning that the number of particles in the simulation is different from numBonds
            std::cout << "Error. " << numBonds << " bonds desired in simulation but " << countBonds << " bonds in " << this->bondFileName;
        }
        
        inputfile.close();
    }
}

void Simulation::readAngleFile()
{
    ifstream inputfile (angleFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numAngles = 0;
        int countAngles = 0;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numAngles;
            }
            else if(linecount > 1)
            {
                int particleI, particleJ, particleK, angleType;
                iss >> particleI >> particleJ >> particleK >> angleType;
                
                if(particleI < 0 || particleI > pinfo.getNumParticles() || particleJ < 0 || particleJ > pinfo.getNumParticles() || particleK < 0 || particleK > pinfo.getNumParticles())
                {
                    // warning about particleI and/or particleJ being out of bounds
                    std::cout << "One of " << particleI << " "  << particleJ << " " << particleK << " is out of bounds, check " << this->angleFileName;
                }
                
                pinfo.addAngle(Angle {particleI, particleJ, particleK, angleType});
                countAngles++;
                //add angle to ainfo
                //int angleTag ainfo.addAngle(particleI, particleJ, particleK);
                
                //add bond to both particleI and particleJ
                //pinfo.addAngle(particleI, angleType);
                //pinfo.addAngle(particleJ, angleType);
                //pinfo.addAngle(particleK, angleType);
            }
            
            linecount++;
        }
        
        if(countAngles != numAngles)
        {
            // warning that the number of angles in the simulation is different from numAngles
            throw std::runtime_error("Number of angles desired in simulation is " + std::to_string(numAngles) + " but " + std::to_string(countAngles) + " angles in " + angleFileName);
        }
        
        inputfile.close();
    }
}

void Simulation::moveParticles()
{
    int numParticles = pinfo.getNumParticles();
    
    #pragma omp parallel for
    for(int i = 0; i < numParticles; i++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(i);
        Coordinate currForce = pinfo.getForceByIndex(i);
        
        double vx, vy, vz;
        
        if(pinfo.getViscosityTypeByIndex(i) == 0)
        {
            // not in nascent fa
            vx = currForce.x / zeta;
            vy = currForce.y / zeta;
            vz = currForce.z / zeta;
        }
        else
        {
            // in nascent fa
            vx = currForce.x / zetaInNascentFA;
            vy = currForce.y / zetaInNascentFA;
            vz = currForce.z / zetaInNascentFA;
        }
        
        /*if(currPosition.z <= FABeadPos.z + FABeadRadius && currPosition.z >= FABeadPos.z - FABeadRadius && currPosition.y <= 0.0)
        {
            Coordinate currDisp = simBox.calcDisplacement(currPosition, FABeadPos);
            currDisp.y = 0.0;
            
            double currDist = currDisp.getMagnitude();
            if(currDist <= FABeadRadius)
            {
                // in nascent fa
                vx = currForce.x / zetaInNascentFA;
                vy = currForce.y / zetaInNascentFA;
                vz = currForce.z / zetaInNascentFA;
            }
            else
            {
                // not in nascent fa
                vx = currForce.x / zeta;
                vy = currForce.y / zeta;
                vz = currForce.z / zeta;
            }
        }
        else
        {
            // not in nascent fa
            vx = currForce.x / zeta;
            vy = currForce.y / zeta;
            vz = currForce.z / zeta;
        }*/
    
                
        Coordinate nextPosition;
        nextPosition.x = currPosition.x + vx * dt;
        nextPosition.y = currPosition.y + vy * dt;
        nextPosition.z = currPosition.z + vz * dt;
        
        if(std::isnan(nextPosition.x))
            throw std::runtime_error("Particle " + std::to_string(i) + " has a NaN position in the x direction");
        
        nextPosition = simBox.periodicWrap(nextPosition);
        
        pinfo.setPosByIndex(i, nextPosition);
	}
}

void Simulation::putAllParticlesInGrid()
{
    particleGrid.clearGrid();
    
    int numParticles = pinfo.getNumParticles();
    for(int i = numParticles-1; i > -1; i--)
    {
        int tempTag = pinfo.getTagAtIndex(i);
        
        // if this tag exists
        if(tempTag >= 0)
        {
            // can put particles in bins twice!!!!
            int newBin = particleGrid.putInGrid(tempTag, pinfo.getPosByIndex(i));
            
            // if this particle is outside the grid remove it here
            // dont try to remove this from the grid since grid is empty right now
            /*if(newBin == -1)
            {
               // cout << "putAllParticlesInGrid remove: " << i << " " << tempTag << " " << newBin << endl;
                //Coordinate currPos = pinfo.getPosByIndex(i);
                
                //cout << "In putAllParticlesInGrid index: " << i << endl;
                removeParticleByIndex(i);
            }*/
        }
    }
    
    // set new particle neighbors
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        pinfo.setNeighborsByIndex(p, particleGrid.getNeighboringTags(pinfo.getTagAtIndex(p)));
        
        //std::vector <int> tempNeighbors = pinfo.getNeighbors(pinfo.getTagAtIndex(p));
    }
    
    std::vector <int> possibleBondNeighborTypes = {0};
    // now set bond neighbors but only for bonds with type 0 (filament bonds are only included for excluded volume)
    if(excludedVolOn)
    {
        //std::vector <int> bondNeighborDesiredOrder;
        
        #pragma omp parallel for
        for(int b = 0; b < pinfo.getNumBonds(); b++)
        {
            pinfo.resetBondNeighborTagsByIndex(b, possibleBondNeighborTypes);
            
            // calc bond distance here between b and neighbors and if they are greater than some value then remove them
            
            Bond bondB = pinfo.getBondByIndex(b);
            Coordinate b_i = pinfo.getPos(bondB.i);
            Coordinate b_j = pinfo.getPos(bondB.j);
                
            std::vector <int> newNeighboringTags = pinfo.getBondNeighborTagsByIndex(b);
            //cout << newNeighboringTags.size() << ", ";
            
            for(int neigh = newNeighboringTags.size()-1; neigh > -1; neigh--)
            {
                Bond neighborBond = pinfo.getBond(newNeighboringTags[neigh]);
                Coordinate neigh_i = pinfo.getPos(neighborBond.i);
                Coordinate neigh_j = pinfo.getPos(neighborBond.j);
                
                double currMinDist = calcMinDistanceCylinders(b_i, b_j, neigh_i, neigh_j);
                
                if(currMinDist > 0.05)
                {
                    // remove this bond from newNeighboringTags
                    newNeighboringTags[neigh] = newNeighboringTags.back();
                    newNeighboringTags.pop_back();
                }
            }
            
           // cout << newNeighboringTags.size() << endl;
            pinfo.setBondNeighborTagsByIndex(b, newNeighboringTags);
            
            /*
            Bond masterBond = pinfo.getBondByIndex(b);
            if(masterBond.bondTypeIndex == 0)
            {
                
                
                /*int bondTag = pinfo.getBondTagAtPos(b);
                
                int tagI = masterBond.i;
                int tagJ = masterBond.j;
                
                std::vector <int> iNeighbors = pinfo.getNeighbors(tagI);
                std::vector <int> jNeighbors = pinfo.getNeighbors(tagJ);
                
                //https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
                //allNeighbors is a set which contains all beads which are either a neighbor of tagJ or tagI
                std::set <int> allNeighbors;
                for(int m = 0; m < iNeighbors.size(); m++)
                    allNeighbors.insert(iNeighbors[m]);
                for(int m = 0; m < jNeighbors.size(); m++)
                    allNeighbors.insert(jNeighbors[m]);
                
                // now loop through allNeighbors to find bond tags associated with these neighbors
                std::set <int> bondNeighborTags;
                
                for(std::set<int>::iterator it = allNeighbors.begin(); it != allNeighbors.end(); ++it)
                {
                    int currNeighborTagI = *it;
                
                    //loop over all bonds connecting to particle currNeighborTagI
                    std::vector <struct HalfBond> secondBondList = pinfo.getBondTags(currNeighborTagI);
                    
                    for(int loopVar = 0; loopVar < secondBondList.size(); loopVar++)
                    {
                        HalfBond tempHalfBond = secondBondList[loopVar];
                        int currNeighborTagJ = tempHalfBond.j;
                        int secondBondTag = tempHalfBond.bondTag;
                        
                        int secondBondType = pinfo.getBond(secondBondTag).bondTypeIndex;
                        
                        bool isNeighboringBond = (tagI == currNeighborTagI) || (tagI == currNeighborTagJ) || (tagJ == currNeighborTagI) || (tagJ == currNeighborTagJ);
                        
                        // only add bonds as neighbors that are type 0 (filament bonds), are of a higher tag (don't want to double count when calculating excluded volume force but this is needed when removing bonds!!)
                        // and are not neighboring
                        if(secondBondType == 0 && !isNeighboringBond)
                        {
                            bondNeighborTags.insert(secondBondTag);
                            
                            
                        }
                    }
                }
                
                std::vector <int> tempVector;
                tempVector.assign(bondNeighborTags.begin(), bondNeighborTags.end());
                
                // dont need this critical statement?
                //#pragma omp critical
                {
                    pinfo.setBondNeighborTagsByIndex(b, tempVector);
                }*/
           // }
        }
    }
    
   // exit(0);
}

void Simulation::calcFilamentForces(int f)
{
    std::deque <int> orgTagList = filaments[f].getAllTags();
    
    int numTags = orgTagList.size();
    int removedTags = 0;
    
    Coordinate ab, bc;
    Coordinate abUnit, bcUnit;
    int aTag, bTag, cTag;
    Coordinate posBTag, posCTag;
    double abMag, bcMag;
    
    for(int p = 0; p < numTags-1; p++)
    {
        // shift by one bead to the left since last loop calculated these parameters
        ab = bc;
        aTag = bTag;
        abUnit = bcUnit;
        abMag = bcMag;
        
        if(p == 0)
        {
            bTag = orgTagList[p];
            posBTag = pinfo.getPos(bTag);
        }
        else
        {
            bTag = cTag;
            posBTag = posCTag;
        }

        cTag = orgTagList[p+1];
        
        posCTag = pinfo.getPos(cTag);
        
        bc = simBox.calcDisplacement(posBTag, posCTag);
        bcMag = bc.getMagnitude();
        bcUnit = bc / bcMag;
        
        // assumes that all bonds in a filament are the same and 
        // the parameters for this bond are located in the first position of vector bondTypes
        double bondForce = bondTypes[0].calcForce(bcMag);
        /*if(fabs(bondForce) > 5.0)
        {
            cout << bondForce << " " << endl;
        }*/
        
        bool brokeDueToBending = false;
        if(p > 0)
        {   
            // calc bending forces
            // dont want to do this if just deleted the previous bond
            Coordinate cb = -bc;
            
            double averageMag = (abMag + bcMag) * 0.5;
            double invAverageMag = 1.0 / averageMag;
            
            double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
            
            Coordinate forceA = -filamentK * invAverageMag * (cb / bcMag / abMag + dotProd / bcMag * (-ab / (abMag*abMag*abMag)));
            Coordinate forceC = -filamentK * invAverageMag * (ab / abMag / bcMag + dotProd / abMag * (-cb / (bcMag*bcMag*bcMag)));
            
            Coordinate forceB = forceA + forceC;
            
            // if forceB mag is greater than some value break filament
            if(false && forceB.getMagnitude() > criticalBendingForceMag)
            {
                //#pragma omp critical
                {
                    // search for bond tag connecting b & c
                    int tempBondTag = -1;
                    std::vector <HalfBond> tempBondList = pinfo.getBondTags(bTag);
                    for(int s = 0; s < tempBondList.size(); s++)
                    {
                        if(tempBondList[s].j == cTag)
                        {
                            tempBondTag = tempBondList[s].bondTag;
                            break;
                        }
                    }
                    
                    if(tempBondTag >= 0)
                    {
                        // flag this bond for removal by removeFlaggedBonds()
                        cout << "Flagging bond for bending removal with tag: " << tempBondTag << ", between particles: " << bTag << " " << cTag << endl;
                        //pinfo.setBondType(tempBondTag, -2);
                        //instead add this tag to the removalList
                        
                        #pragma omp critical
                        {
                            removalList.push_back(tempBondTag);
                        }
                    }
                    else
                    {            
                        throw std::runtime_error("Could not find bondTag between particles with tags: " + std::to_string(bTag) + " " + std::to_string(cTag));
                    }
                    brokeDueToBending = true;
                }
            }
            else
            {
                pinfo.addForce(aTag, forceA);
                pinfo.addForce(bTag, -forceB);
                pinfo.addForce(cTag, forceC);
                
                Coordinate tempStress;
                
                tempStress.x = 1.0/3.0 * ((ab.x+posBTag.x)*forceA.x - posBTag.x*forceB.x + (cb.x+posBTag.x)*forceC.x);
                tempStress.y = 1.0/3.0 * ((ab.y+posBTag.y)*forceA.y - posBTag.y*forceB.y + (cb.y+posBTag.y)*forceC.y);
                tempStress.z = 1.0/3.0 * ((ab.z+posBTag.z)*forceA.z - posBTag.z*forceB.z + (cb.z+posBTag.z)*forceC.z);
                
                pinfo.addStress(aTag, tempStress, 1);
                pinfo.addStress(bTag, tempStress, 1); 
                pinfo.addStress(cTag, tempStress, 1);
            }
        }
        
        Coordinate tempForce = bondForce * bcUnit;
        Coordinate tempStress;
            
        tempStress.x = 0.5 * (tempForce.x * bc.x);
        tempStress.y = 0.5 * (tempForce.y * bc.y);
        tempStress.z = 0.5 * (tempForce.z * bc.z);
        
        double tempStressMag = sqrt(tempStress.x*tempStress.x + tempStress.y*tempStress.y + tempStress.z*tempStress.z);
        
        // this search seems to add 34 seconds per 10 simulation seconds (387.393->422.066)
        // break filament into two filaments if bondForce is too large?
        // but cant do this here since looping over filaments?
        // search for bond tag connecting b & c
        int tempBondTag = -1;
        std::vector <HalfBond> tempBondList = pinfo.getBondTags(bTag);
        for(int s = 0; s < tempBondList.size(); s++)
        {
            if(tempBondList[s].j == cTag)
            {
                tempBondTag = tempBondList[s].bondTag;
                break;
            }
        }
        
        double currBondDestructionTime = 0.0;
        if(tempBondTag >= 0)
        {
            currBondDestructionTime = pinfo.getBond(tempBondTag).destructionTime;
        }
        else
        {            
            throw std::runtime_error("Could not find bondTag between particles with tags: " + std::to_string(bTag) + " " + std::to_string(cTag));
        }
        

        if(currSimTime >= currBondDestructionTime || false && bondForce < criticalBondForceMag && brokeDueToBending == false)
        {
            //cout << "Flagging bond for removal with tag: " << tempBondTag << ", between particles: " << bTag << " " << cTag << endl;
            //pinfo.setBondType(tempBondTag, -1);
            //instead add this tag to the removalList
            
            #pragma omp critical
            {
                removalList.push_back(tempBondTag);
            }
        }
        else if(brokeDueToBending == false)
        {            
            //Coordinate tempForce = bondForce * bcUnit;
            pinfo.addForce(bTag, tempForce);
            pinfo.addForce(cTag, -tempForce);
            
            // stress calculations
            /*Coordinate tempStress;
            
            tempStress.x = 0.5 * (tempForce.x * bc.x);
            tempStress.y = 0.5 * (tempForce.y * bc.y);
            tempStress.z = 0.5 * (tempForce.z * bc.z);*/
            
            pinfo.addStress(bTag, tempStress, 0);
            pinfo.addStress(cTag, tempStress, 0);
        }
    }    
}

void Simulation::calcBondForce(Bond b)
{
    Coordinate bondForceVector = getBondForce(b);
    
    pinfo.addForce(b.i, bondForceVector);
    pinfo.addForce(b.j, -bondForceVector);
}

Coordinate Simulation::getBondForce(Bond b)
{
    Coordinate i = pinfo.getPos(b.i);
    Coordinate j = pinfo.getPos(b.j);
    
    int bondType = b.bondTypeIndex;
    
    Coordinate ij = simBox.calcDisplacement(i, j);
    
    return getBondForce(ij, bondType);
}

// calculates bond force for a displacement of ij for a bond of type bondType
Coordinate Simulation::getBondForce(Coordinate ij, int bondType)
{
    double ijMag = ij.getMagnitude();
    if(ijMag == 0.0)
    {
        return Coordinate {0.0, 0.0, 0.0};
    }
    
    Coordinate ijUnit = ij / ijMag;
    
    // assumes that all bonds in a filament are the same and 
    // the parameters for this bond are located in the first position of vector bondTypes
    double bondForce = bondTypes[bondType].calcForce(ijMag);
    
    return bondForce * ijUnit;
}

void Simulation::calcFilamentAngleForce(Angle ang)
{
    Coordinate a = pinfo.getPos(ang.i);
    Coordinate b = pinfo.getPos(ang.j);
    Coordinate c = pinfo.getPos(ang.k);
    
    Coordinate ab = simBox.calcDisplacement(a, b);
    double abMag = ab.getMagnitude();
    Coordinate abUnit = ab / abMag;
    
    Coordinate cb = simBox.calcDisplacement(c, b);
    double cbMag = cb.getMagnitude();
    Coordinate cbUnit = cb / cbMag;
    
    
    double averageMag = (abMag + cbMag) * 0.5;
    double invAverageMag = 1.0 / averageMag;
    
    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
    
    Coordinate forceA = -filamentK * invAverageMag * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
    Coordinate forceC = -filamentK * invAverageMag * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));
    
    Coordinate forceB = forceA + forceC;
    
    pinfo.addForce(ang.i, forceA);
    pinfo.addForce(ang.j, -forceB);
    pinfo.addForce(ang.k, forceC);
}

void Simulation::calcAngleForce(Angle ang)
{
    double theta0 = cos(1.22173);
    double angleK = branchAngleK;
    
    Coordinate a = pinfo.getPos(ang.i);
    Coordinate b = pinfo.getPos(ang.j);
    Coordinate c = pinfo.getPos(ang.k);
    
    Coordinate ab = simBox.calcDisplacement(a, b);
    double abMag = ab.getMagnitude();
    Coordinate abUnit = ab / abMag;
    
    Coordinate cb = simBox.calcDisplacement(c, b);
    double cbMag = cb.getMagnitude();
    Coordinate cbUnit = cb / cbMag;
    
    double tempDot = abUnit.x*cbUnit.x + abUnit.y*cbUnit.y + abUnit.z*cbUnit.z;
    
    if(tempDot < -1.0)
    {
        tempDot = -1.0;
    }
    else if(tempDot > 1.0)
    {
        tempDot = 1.0;
    }
    
    double theta = acos(tempDot);
    
    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
    
    Coordinate forceA = -2.0*angleK * (tempDot - theta0) * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
    Coordinate forceC = -2.0*angleK * (tempDot - theta0) * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));
    
    Coordinate forceB = forceA + forceC;

    pinfo.addForce(ang.i, forceA);
    pinfo.addForce(ang.j, -forceB);
    pinfo.addForce(ang.k, forceC);
}

// returns forces in the input parameters
Coordinate Simulation::calcExcludedVolumeForces(Coordinate& ri, Coordinate& rj, Coordinate& ri2, Coordinate& rj2)
{
    Coordinate tempStress = {0.0, 0.0, 0.0};
    
    /*ri = Coordinate {0.0, 0.0, 0.0};
    rj = Coordinate {0.0, 0.0, 0.0};
    ri2 = Coordinate {0.0, 0.0, 0.0};
    rj2 = Coordinate {0.0, 0.0, 0.0};
        
    return tempStress;*/
    
    Coordinate R_k = simBox.calcDisplacement(ri, rj);
    
    Coordinate P_k = 0.5 * R_k + rj;
    P_k = simBox.periodicWrap(P_k);
    
    double R_kSquared = R_k*R_k;

    Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

    Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);
    
    double R_lSquared = R_l*R_l;
    
    double R_klSquared = R_k*R_l;
    
    Coordinate firstTerm = simBox.calcDisplacement(P_k, P_l);
    
    Coordinate secondTerm = simBox.periodicWrap(R_lSquared * R_k - R_klSquared*R_l);
    
    double tk = (firstTerm) * (secondTerm);
    tk = tk / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tk > 0.5)
        tk = 0.5;
    else if(tk < -0.5)
        tk = -0.5;
    
    // does this work?
    firstTerm = -firstTerm;
    
    secondTerm = simBox.periodicWrap(R_kSquared * R_l - R_klSquared*R_k);
    
    double tl = (firstTerm) * (secondTerm);
    tl = tl / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tl > 0.5)
        tl = 0.5;
    else if(tl < -0.5)
        tl = -0.5;
    
    Coordinate closestPointOnK = simBox.periodicWrap(P_k + tk * R_k);
    
    Coordinate closestPointOnL = simBox.periodicWrap(P_l + tl * R_l);

    // points towards k from l rather than the reverse
    Coordinate d_kl = simBox.calcDisplacement(closestPointOnK, closestPointOnL);
    
    double minDist = d_kl.getMagnitude();					

    // calc potential regardless of minDist?
    if(minDist < rRepulsiveInteraction)
    {
        // kr is too large by ~2 orders of magnitude compared to Screiber, Stewart and Duke 2010
        double forceMag = 0.0;
        
        double kr = -100.0;
        // excluded vol
        forceMag = kr * (minDist - rRepulsiveInteraction);

        Coordinate Fr = forceMag * d_kl.getUnitCoord();
        
        double length_k = R_k.getMagnitude();
        
        double firstSegmentLength = (tk + 0.5) * length_k;
        double secondSegmentLength = length_k - firstSegmentLength;
        
        Coordinate Falpha = firstSegmentLength / length_k * (Fr);
        Coordinate Fbeta = secondSegmentLength / length_k * (Fr);
        
        double length_l = R_l.getMagnitude();
        
        firstSegmentLength = (tl + 0.5) * length_l;
        secondSegmentLength = length_l - firstSegmentLength;
        
        Coordinate Falpha2 = firstSegmentLength / length_l * (-Fr);
        Coordinate Fbeta2 = secondSegmentLength / length_l * (-Fr);
        
        ri = Falpha;
        rj = Fbeta;
        ri2 = Falpha2;
        rj2 = Fbeta2;
                                         
        Coordinate R3 = simBox.calcDisplacement(ri2, rj);
        
        Coordinate R4 = simBox.calcDisplacement(ri2, rj);
        
        tempStress.x = 0.25 * ((R_k.x+rj.x)*Falpha.x + rj.x*Fbeta.x + (R3.x+rj.x)*Falpha2.x + (R4.x+rj.x)*Fbeta2.x);
        tempStress.y = 0.25 * ((R_k.y+rj.y)*Falpha.y + rj.y*Fbeta.y + (R3.y+rj.y)*Falpha2.y + (R4.x+rj.x)*Fbeta2.x);
        tempStress.z = 0.25 * ((R_k.z+rj.z)*Falpha.z + rj.z*Fbeta.z + (R3.z+rj.z)*Falpha2.z + (R4.x+rj.x)*Fbeta2.x);
    }
    else
    {
        ri = Coordinate {0.0, 0.0, 0.0};
        rj = Coordinate {0.0, 0.0, 0.0};
        ri2 = Coordinate {0.0, 0.0, 0.0};
        rj2 = Coordinate {0.0, 0.0, 0.0};
    }
    
    return tempStress;
}


double Simulation::calcMinDistanceCylinders(Coordinate ri, Coordinate rj, Coordinate ri2, Coordinate rj2)
{
    Coordinate R_k = simBox.calcDisplacement(ri, rj);
    
    Coordinate P_k = 0.5 * R_k + rj;
    P_k = simBox.periodicWrap(P_k);
    
    double R_kSquared = R_k*R_k;

    Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

    Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);
    
    double R_lSquared = R_l*R_l;
    
    double R_klSquared = R_k*R_l;
    
    Coordinate firstTerm = simBox.calcDisplacement(P_k, P_l);
    
    Coordinate secondTerm = simBox.periodicWrap(R_lSquared * R_k - R_klSquared*R_l);
    
    double tk = (firstTerm) * (secondTerm);
    tk = tk / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tk > 0.5)
        tk = 0.5;
    else if(tk < -0.5)
        tk = -0.5;
    
    // does this work?
    firstTerm = -firstTerm;
    
    secondTerm = simBox.periodicWrap(R_kSquared * R_l - R_klSquared*R_k);
    
    double tl = (firstTerm) * (secondTerm);
    tl = tl / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tl > 0.5)
        tl = 0.5;
    else if(tl < -0.5)
        tl = -0.5;
    
    Coordinate closestPointOnK = simBox.periodicWrap(P_k + tk * R_k);
    
    Coordinate closestPointOnL = simBox.periodicWrap(P_l + tl * R_l);

    // points towards k from l rather than the reverse
    Coordinate d_kl = simBox.calcDisplacement(closestPointOnK, closestPointOnL);
    
    double minDist = d_kl.getMagnitude();					
    
    return minDist;
}

double Simulation::calcForces()
{
    int numParticles = pinfo.getNumParticles();

    // make tempForces vectors larger if needed and set new elements to zero
    // tempStress currently only holds one type of stress (bond)!
    if(tempForces[0].size() != numParticles)
    {
        #pragma omp parallel for
        for(int thread = 0; thread < omp_get_max_threads(); thread++)
        {
            int oldSize = tempForces[thread].size();
            
            tempForces[thread].resize(pinfo.getNumParticles());
            tempStresses[thread].resize(pinfo.getNumParticles());
            
            for(int i = oldSize; i < pinfo.getNumParticles(); i++)
            {
                tempForces[thread][i] = Coordinate {0.0, 0.0, 0.0};
                tempStresses[thread][i] = Coordinate {0.0, 0.0, 0.0};
            }
        }
    }
    
	// calculate other forces using a scheduler based on the bins they are in
    // or by using tempForces vector (this uses more memory)
    #pragma omp parallel for
    for(int p = 0; p < numParticles; p++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(p);
        
        int currTag = pinfo.getTagAtIndex(p);
        std::vector <int> iNeighbors = pinfo.getNeighbors(currTag);
        
        // loop over extra bonded particles
        std::vector <struct HalfBond> tempBondList = pinfo.getBondTagsByIndex(p);
        for(int b = 0; b < tempBondList.size(); b++)
        {
            Bond firstBond = pinfo.getBond(tempBondList[b].bondTag);
            
            // this was wrong in ASCB version
            if(currTag < tempBondList[b].j)
            {
                // bondTag is  being used here as bondType and bond list is not used at all in pinfo
                //calcBondForce(Bond {currTag, tempBondList[b].j, tempBondList[b].bondTag});
                //Bond tempBond = {currTag, tempBondList[b].j, 1};
                
                Coordinate ri = pinfo.getPos(currTag);
                Coordinate rj = pinfo.getPos(firstBond.j);
                Coordinate rij = simBox.calcDisplacement(ri, rj);
                
                // if this bond is a crosslink bond then calculate bond forces on it (but no excluded forces)
                if(firstBond.bondTypeIndex > 0)
                {
                    double rijMag = rij.getMagnitude();
                    
                    //double bondLifetime = currSimTime - firstBond.creationTime;
                    double bondDestructionTime = firstBond.destructionTime;
                    
                    /*#pragma omp critical
                    {
                        cout << bondDestructionTime << " " << currSimTime << " " << firstBond.bondTypeIndex << endl;
                    }*/
                    
                    // temporary bonds (type 2) have lifetimes of 1.0
                    //if(firstBond.bondTypeIndex == 1 || (firstBond.bondTypeIndex == 2 && bondDestructionTime < currSimTime) || (severCrosslinksWithExtension && rijMag < 2.0*maxCrosslinkBondDist))
                    
                    if((currSimTime > bondDestructionTime) || (severCrosslinksWithExtension && rijMag > 2.0*maxCrosslinkBondDist))
                    {
                        // flag this bond for removal by removeFlaggedBonds()
                        //pinfo.setBondType(tempBondList[b].bondTag, -1);
                        //instead add this tag to the removalList
                        
                        #pragma omp critical
                        {
                            removalList.push_back(tempBondList[b].bondTag);
                        }
                    }
                    else
                    {
                        Coordinate rijUnit = rij / rijMag;
                        
                        Coordinate bondForceVector = bondTypes[firstBond.bondTypeIndex].calcForce(rijMag) * rijUnit;
                    
                        //int currPosI = pinfo.getIndexOfTag(currTag);
                        int currPosJ = pinfo.getIndexOfTag(tempBondList[b].j);
                    
                        tempForces[omp_get_thread_num()][p] = tempForces[omp_get_thread_num()][p] + bondForceVector;
                        tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] - bondForceVector;
                        
                        //stress calculations
                        Coordinate tempStress;
                        
                        tempStress.x = 0.5 * (bondForceVector.x * rij.x);
                        tempStress.y = 0.5 * (bondForceVector.y * rij.y);
                        tempStress.z = 0.5 * (bondForceVector.z * rij.z);
                        
                        tempStresses[omp_get_thread_num()][p] = tempStresses[omp_get_thread_num()][p] + tempStress;
                        tempStresses[omp_get_thread_num()][currPosJ] = tempStresses[omp_get_thread_num()][currPosJ] + tempStress;
                    }
                }
                
            }
        }
        
        /*
        // loop over extra angle particles (these might be too far away?)
        std::vector <struct Angle> tempAngleList = pinfo.getAngleTagsByIndex(i);
        for(int a = 0; a < tempAngleList.size(); a++)
        {
            if(currTag == tempAngleList[a].i)
            {
                if(tempAngleList[a].angleTypeIndex == 0)
                {
                    calcFilamentAngleForce(tempAngleList[a]);
                }
                else
                {
                    calcAngleForce(tempAngleList[a]);
                }
            }
        }*/
        
	}
    
    if(excludedVolOn)
    {
        // remove bonds before excluded volume calculations??? or in removeFlaggedBonds?
        int currNumBonds = pinfo.getNumBonds();
        
        std::vector <int> countPerThread (omp_get_max_threads());
        //for(int i = 0; i < countPerThread.size(); i++)
        //    countPerThread[i] = 0;
        
        
        // bonds of lower indexes have larger numbers of neighbors on average so the work per thread will not be balanced
        #pragma omp parallel for schedule(static, 1)
        for(int b = 0; b < currNumBonds; b++)
        {            
            std::vector <int> currNeighborBondTags = pinfo.getBondNeighborTagsByIndex(b);
            Bond firstBond = pinfo.getBondByIndex(b);
            
            int currPosI, currPosJ;
            Coordinate ri, rj;
            
            try
            {
                currPosI = pinfo.getIndexOfTag(firstBond.i);
                currPosJ = pinfo.getIndexOfTag(firstBond.j);
                
                ri = pinfo.getPosByIndex(currPosI);
                rj = pinfo.getPosByIndex(currPosJ);
            }
            catch (const std::runtime_error& error)
            {
                #pragma omp critical
                {
                    cout << b << " " << firstBond.i << " " << firstBond.j << " " << currNumBonds << " " << currPosI << " " << currPosJ << endl;
                    exit(0);
                }
            }
                
            // all bond distances must have been calculated before, store these somehow for use here?
            Coordinate R_k = simBox.calcDisplacement(ri, rj);
            double R_kSquared = R_k*R_k;
            
            Coordinate P_k = 0.5 * R_k + rj;
            P_k = simBox.periodicWrap(P_k);
            
            //int calcCount = 0;
            //int withinRangeCount = 0;
            
            for(int neigh = 0; neigh < currNeighborBondTags.size(); neigh++)
            {
                // if tag of neighboring bond is less than tag at pos b then calculate force
                if(currNeighborBondTags[neigh] < pinfo.getBondTagAtPos(b))
                {
                    //calcCount++;
                    
                   // countPerThread[omp_get_thread_num()] += 1;
                    
                    /*Bond secondBond = pinfo.getBond(currNeighborBondTags[neigh]);
                    
                    Coordinate ri = pinfo.getPosByIndex(currPosI);
                    Coordinate rj = pinfo.getPosByIndex(currPosJ);
                    
                    int currPosI2 = pinfo.getIndexOfTag(secondBond.i);
                    int currPosJ2 = pinfo.getIndexOfTag(secondBond.j);
                    
                    Coordinate ri2 = pinfo.getPosByIndex(currPosI2);
                    Coordinate rj2 = pinfo.getPosByIndex(currPosJ2);
                    
                    Coordinate tempStress = calcExcludedVolumeForces(ri, rj, ri2, rj2);

                    tempForces[omp_get_thread_num()][currPosI] = tempForces[omp_get_thread_num()][currPosI] + ri;
                    tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] + rj;
                    tempForces[omp_get_thread_num()][currPosI2] = tempForces[omp_get_thread_num()][currPosI2] + ri2;
                    tempForces[omp_get_thread_num()][currPosJ2] = tempForces[omp_get_thread_num()][currPosJ2] + rj2;*/
                    
                    
                    Coordinate tempStress = {0.0, 0.0, 0.0};
                    
                    Bond secondBond = pinfo.getBond(currNeighborBondTags[neigh]);
                    
                    int currPosI2, currPosJ2;
                    Coordinate ri2, rj2;
            
                    try
                    {
                        currPosI2 = pinfo.getIndexOfTag(secondBond.i);
                        currPosJ2 = pinfo.getIndexOfTag(secondBond.j);
                        
                        ri2 = pinfo.getPosByIndex(currPosI2);
                        rj2 = pinfo.getPosByIndex(currPosJ2);
                    }
                    catch(const runtime_error& error)
                    {
                        #pragma omp critical
                        {
                            cout << "second" << " " << secondBond.i << " " << secondBond.j << " " << currPosI2 << " " << currPosJ2 << endl;
                            exit(0);
                        }
                    }

                    Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

                    Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);
                    
                    double R_lSquared = R_l*R_l;
                    
                    double R_klSquared = R_k*R_l;
                    
                    Coordinate firstTerm = simBox.calcDisplacement(P_k, P_l);
                    
                    Coordinate secondTerm = simBox.periodicWrap(R_lSquared * R_k - R_klSquared*R_l);
                    
                    double tk = (firstTerm) * (secondTerm);
                    tk = tk / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
                    
                    if(tk > 0.5)
                        tk = 0.5;
                    else if(tk < -0.5)
                        tk = -0.5;
                    
                    // does this work?
                    firstTerm = -firstTerm;
                    
                    secondTerm = simBox.periodicWrap(R_kSquared * R_l - R_klSquared*R_k);
                    
                    double tl = (firstTerm) * (secondTerm);
                    tl = tl / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
                    
                    if(tl > 0.5)
                        tl = 0.5;
                    else if(tl < -0.5)
                        tl = -0.5;
                    
                    Coordinate closestPointOnK = simBox.periodicWrap(P_k + tk * R_k);
                    
                    Coordinate closestPointOnL = simBox.periodicWrap(P_l + tl * R_l);

                    // points towards k from l rather than the reverse
                    Coordinate d_kl = simBox.calcDisplacement(closestPointOnK, closestPointOnL);
                    
                    double minDist = d_kl.getMagnitude();					


                    
                    // calc potential regardless of minDist?
                    if(minDist < rRepulsiveInteraction)
                    {
                       // withinRangeCount++;
                        double forceMag = 0.0;
                        
                        double kr = -1690.0;
                        // excluded vol
                        forceMag = kr * (minDist - rRepulsiveInteraction);

                        Coordinate Fr = forceMag * d_kl.getUnitCoord();
                        
                        double length_k = R_k.getMagnitude();
                        
                        double firstSegmentLength = (tk + 0.5) * length_k;
                        double secondSegmentLength = length_k - firstSegmentLength;
                        
                        Coordinate Falpha = firstSegmentLength / length_k * (Fr);
                        Coordinate Fbeta = secondSegmentLength / length_k * (Fr);
                        
                        double length_l = R_l.getMagnitude();
                        
                        firstSegmentLength = (tl + 0.5) * length_l;
                        secondSegmentLength = length_l - firstSegmentLength;
                        
                        Coordinate Falpha2 = firstSegmentLength / length_l * (-Fr);
                        Coordinate Fbeta2 = secondSegmentLength / length_l * (-Fr);
                        
                        if(currPosI >= tempForces[omp_get_thread_num()].size() || currPosJ >= tempForces[omp_get_thread_num()].size() || currPosI2 >= tempForces[omp_get_thread_num()].size() || currPosJ2 >= tempForces[omp_get_thread_num()].size())
                        {
                            cout << "too big for tempForces: " << currPosI << " " << currPosJ << " " << currPosI2 << " " << currPosJ2 << " " << pinfo.getNumParticles() << " " << tempForces[omp_get_thread_num()].size() << endl;
                            exit(0);
                        }
                        
                        if(currPosI < 0 || currPosJ < 0 || currPosI2 < 0 || currPosJ2 < 0)
                        {
                            cout << "out of bounds of tempForces (-): " << currPosI << " " << currPosJ << " " << currPosI2 << " " << currPosJ2 << " " << pinfo.getNumParticles() << " " << tempForces[omp_get_thread_num()].size() << endl;
                            exit(0);
                        }
                        
                        tempForces[omp_get_thread_num()][currPosI] = tempForces[omp_get_thread_num()][currPosI] + Falpha;
                        tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] + Fbeta;
                        tempForces[omp_get_thread_num()][currPosI2] = tempForces[omp_get_thread_num()][currPosI2] + Falpha2;
                        tempForces[omp_get_thread_num()][currPosJ2] = tempForces[omp_get_thread_num()][currPosJ2] + Fbeta2;
                        
                        
                                                         
                        /*Coordinate R3 = simBox.calcDisplacement(ri2, rj);
                        
                        Coordinate R4 = simBox.calcDisplacement(ri2, rj);
                        
                        tempStress.x = 1.0/4.0 * ((R_k.x+rj.x)*Falpha.x + rj.x*Fbeta.x + (R3.x+rj.x)*Falpha2.x + (R4.x+rj.x)*Fbeta2.x);
                        tempStress.y = 1.0/4.0 * ((R_k.y+rj.y)*Falpha.y + rj.y*Fbeta.y + (R3.y+rj.y)*Falpha2.y + (R4.x+rj.x)*Fbeta2.x);
                        tempStress.z = 1.0/4.0 * ((R_k.z+rj.z)*Falpha.z + rj.z*Fbeta.z + (R3.z+rj.z)*Falpha2.z + (R4.x+rj.x)*Fbeta2.x);
                        
                        tempStresses[omp_get_thread_num()][currPosI] = tempStresses[omp_get_thread_num()][currPosI] + tempStress;
                        tempStresses[omp_get_thread_num()][currPosJ] = tempStresses[omp_get_thread_num()][currPosJ] + tempStress;
                        tempStresses[omp_get_thread_num()][currPosI2] = tempStresses[omp_get_thread_num()][currPosI2] + tempStress;
                        tempStresses[omp_get_thread_num()][currPosJ2] = tempStresses[omp_get_thread_num()][currPosJ2] + tempStress;
                    */
                    }
                }
            }
            
            
            /*#pragma omp critical
            {
                cout << (double)withinRangeCount / (double)calcCount * 100.0 << " " << withinRangeCount << " " << calcCount << endl;
            }*/
        }
    }
    
   /* if(excludedVolOn)
    {
        #pragma omp parallel for
        for(int i = 0; i < excludedVolumeList.size(); i++)
        {
            int tagI = excludedVolumeList[i][0];
            int tagJ = excludedVolumeList[i][1];
            int tagI2 = excludedVolumeList[i][2];
            int tagJ2 = excludedVolumeList[i][3];
            
            Coordinate ri = pinfo.getPos(tagI);
            Coordinate rj = pinfo.getPos(tagJ);
            Coordinate ri2 = pinfo.getPos(tagI2);
            Coordinate rj2 = pinfo.getPos(tagJ2);
            
            Coordinate tempStress = calcExcludedVolumeForces(ri, rj, ri2, rj2);
            
            int currPosI = pinfo.getIndexOfTag(tagI);
            int currPosJ = pinfo.getIndexOfTag(tagJ);
            
            int currPosI2 = pinfo.getIndexOfTag(tagI2);
            int currPosJ2 = pinfo.getIndexOfTag(tagJ2);
            
            tempForces[omp_get_thread_num()][currPosI] = tempForces[omp_get_thread_num()][currPosI] + ri;
            tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] + rj;
            tempForces[omp_get_thread_num()][currPosI2] = tempForces[omp_get_thread_num()][currPosI2] + ri2;
            tempForces[omp_get_thread_num()][currPosJ2] = tempForces[omp_get_thread_num()][currPosJ2] + rj2;
        }
    }*/
            
            
    
    double simBoxYMax = 0.1;
    
    
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        // calculate thermal forces
        std::mt19937& threadEngine = generators[omp_get_thread_num()];
        
        Coordinate currPos = pinfo.getPosByIndex(p);

        bool inNascentFA = false;
        
        if(currPos.z <= FABeadPos.z + FABeadRadius && currPos.z >= FABeadPos.z - FABeadRadius - FABeadLength && currPos.y <= 0.0)
        {
            if(currPos.z >= FABeadPos.z)
            {
                // in top part of long nascent FA
                Coordinate currDisp = simBox.calcDisplacement(currPos, FABeadPos);
                currDisp.y = 0.0;
                
                double currDist = currDisp.getMagnitude();
                if(currDist <= FABeadRadius)
                {
                    // in nascent fa
                    //Coordinate c = {thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine)};
                    
                    inNascentFA = true;
                }
            }
            else if (currPos.z >= FABeadPos.z - FABeadLength)
            {
                // in middle part of nascent FA
                if(fabs(currPos.x) <= FABeadRadius)
                {
                    // in nascent FA
                    inNascentFA = true;
                }
            }
            else
            {
                // in bottom part of nascent FA
                Coordinate bottomPos = FABeadPos;
                bottomPos.z = bottomPos.z - FABeadLength;
                
                Coordinate currDisp = simBox.calcDisplacement(currPos, bottomPos);
                currDisp.y = 0.0;
                
                double currDist = currDisp.getMagnitude();
                if(currDist <= FABeadRadius)
                {
                    inNascentFA = true;
                }
            }
        }
        
        if(inNascentFA)
        {
            if(pinfo.getViscosityTypeByIndex(p) == 0)
            {
                // not listed as being in nascent focal adhesion region, see if it binds
                double randomVal = dis(threadEngine);
                
                if(randomVal < nascentFABindProbability)
                {
                    pinfo.setViscosityTypeByIndex(p,1);
                }
            }
            else
            {
                // listed as being in nascent focal adhesion region, see if it unbinds
                double randomVal = dis(threadEngine);
                
                if(randomVal < nascentFAUnbindProbability)
                {
                    pinfo.setViscosityTypeByIndex(p,0);
                }
            }
        }
        else
        {
            // not in nascent fa
            if(pinfo.getViscosityTypeByIndex(p) != 0)
            {
                // listed as being in nascent focal adhesion region, unbind it
                //double randomVal = dis(threadEngine);
                
                //if(randomVal < nascentFAUnbindProbability)
                {
                    pinfo.setViscosityTypeByIndex(p,0);
                }
            }
        }
            
        if(temperature > 0.0 && thermalForcesOn == true)
        {
            // set thermal force
            if(pinfo.getViscosityTypeByIndex(p) == 0)
            {
                Coordinate c = {thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine)};
                pinfo.setForceByIndex(p, c);
            }
            else
            {
                Coordinate c = {thermal_force_inNascentFA*gaussianDis(threadEngine), thermal_force_inNascentFA*gaussianDis(threadEngine), thermal_force_inNascentFA*gaussianDis(threadEngine)};
                pinfo.setForceByIndex(p, c);
            }
            
            
        }
        else
        {
            pinfo.setForceByIndex(p, Coordinate {0.0, 0.0, 0.0});
        }
        
        // zero stresses for all particles
        pinfo.setStressByIndex(p, Coordinate {0.0, 0.0, 0.0}, 0);
        pinfo.setStressByIndex(p, Coordinate {0.0, 0.0, 0.0}, 1);
        
        // boundary force in y-direction
        //currPos = pinfo.getPosByIndex(p);
        if(fabs(currPos.y) >= simBoxYMax)
        {
            pinfo.addForceByIndex(p, 1.0 * Coordinate {0.0, -1.0, 0.0} * copysign(1.0, currPos.y));
        }
        
        if(backPullingForceOn)
        {
            if(currPos.z < -3.75)
            {
                pinfo.addForceByIndex(p, pullForceMag*2 * Coordinate {0.0, 0.0, -1.0});
            }
        }
        
        
        // sum up all entries in tempForces which were split across the threads
        Coordinate totalForce = {0.0, 0.0, 0.0};
        Coordinate totalStress = {0.0, 0.0, 0.0};
        for(int thread = 0; thread < omp_get_max_threads(); thread++)
        {
            Coordinate tempForce = tempForces[thread][p];
            Coordinate tempStress = tempStresses[thread][p];
            
            if(tempForce.x != 0.0 || tempForce.y != 0.0 || tempForce.z != 0.0)
            {
                totalForce = totalForce + tempForce;
                totalStress = totalStress + tempStress;
                
                tempForces[thread][p] = Coordinate {0.0, 0.0, 0.0};
                tempStresses[thread][p] = Coordinate {0.0, 0.0, 0.0};
            }
        }
        pinfo.addForceByIndex(p, totalForce);
        pinfo.addStressByIndex(p, totalStress, 0);
        
        // uniform pulling force
        if(uniformPullingForceOn)
        {
            pinfo.addForceByIndex(p, Coordinate {0.0, 0.0, -pullForceMag});
        }
    }
    
    double totalLeadingEdgeForce = 0.0;
    // calculate filament forces (intra bond and angle forces)
    #pragma omp parallel for reduction (+:totalLeadingEdgeForce)
    for(int f = 0; f < filaments.size(); f++)
    {
        calcFilamentForces(f);
        
        // push whole filament if it is near the leading edge
        totalLeadingEdgeForce += calcLeadingEdgeForce(f);
    }
    
    return totalLeadingEdgeForce;
}

//removes bonds marked for removal, returns the number of crosslinks removed
int Simulation::removeFlaggedBonds()
{
   // int numBondsRemoved = 0;
    
    // will making this parallel work?
    //#pragma omp parallel for
    /*for(int b = pinfo.getNumBonds()-1; b > -1; b--)
    {        
        int currType = pinfo.getBondByIndex(b).bondTypeIndex;
        int currTag = pinfo.getBondTagAtPos(b);
        Bond currBond = pinfo.getBond(currTag);
        if(currType < 0)
        {
            //#pragma omp critical
            {
                if(currType == -1)
                    numBondsRemoved++;
                
                //cout << "removed bond: " << b << " " << currTag << " " << currType << " " << currBond.creationTime << " " << currBond.destructionTime << endl;
                removeBondByIndex(b);
            }
        }
    }
    
    return numBondsRemoved;*/
    
    
    
    int numBondsRemoved = removalList.size();
    
    while(removalList.size() > 0)
    {
        int bTag = removalList.back();
        removalList.pop_back();
        int currIndex = pinfo.getBondPosOfTag(bTag);
        
        if(currIndex >= 0)
        {
            // this bond exists
            //int currType = pinfo.getBondByIndex(currIndex).bondTypeIndex;
            //if(currType < 0)
            {
                //if(currType == -1)
                numBondsRemoved++;
                
                removeBondByIndex(currIndex);
            }
        }
    }
    
    return numBondsRemoved;
}

int Simulation::addTemporaryCrosslinks(int nCrosslinks)
{
     if(nCrosslinks <= 0)
        return 0;
    
    // search for particles between minCrosslinkBondDist and maxCrosslinkBondDist
    // should generate this list when doing excluded volume calculations and using the grid
    std::vector <Bond> candidateBonds;
    // add up to nCrosslinks new crosslinks from z down
    
    #pragma omp parallel for
    for(int i = 0; i < pinfo.getNumParticles(); i++)
    {
        std::mt19937& threadEngine = generators[omp_get_thread_num()];
        
        int particleITag = pinfo.getTagAtIndex(i);
        Coordinate posI = pinfo.getPosByIndex(i);
        
        std::vector <int> tempNeighbors = pinfo.getNeighbors(particleITag);
        
        for(int neighPos = 0; neighPos < tempNeighbors.size(); neighPos++)
        {
            int particleJTag = tempNeighbors[neighPos];
            
            // only want to check pairs between particles once
            if(particleITag < particleJTag)
            {
                Coordinate posJ = pinfo.getPos(particleJTag);
                
                Coordinate uij = simBox.calcDisplacement(posI, posJ);
            
                double distMag = uij.getMagnitude();
                
                if(distMag >= minCrosslinkBondDist && distMag <= maxCrosslinkBondDist)
                {
                    Bond newBond = {particleITag, particleJTag, 2};
                    newBond.creationTime = currSimTime;
                    
                    double randomVal = 0.0;
                    do
                    {
                        randomVal = dis(threadEngine);
                    }
                    while(randomVal == 0.0);
                    
                    newBond.destructionTime = currSimTime -1.0 * log(randomVal);
                    
                    // checkBondExists checks if any bond between particleITag and particleJTag exists
                    bool bondAlreadyExists = pinfo.checkBondExists(newBond);
                    
                    // no actin bond either (just define a new function checkAnyBondExists that disregards type
                    //Bond newBondTwo = {particleITag, particleJTag, 0};
                    //bool bondAlreadyExistsTwo = pinfo.checkBondExists(newBondTwo);
                    
                    if(bondAlreadyExists == false)
                    {
                        double zpos = 0.5*uij.z + posJ.z;
                        
                        #pragma omp critical
                        {
                            //cout << newBond.creationTime << " " << newBond.destructionTime << endl;
                            candidateBonds.push_back(newBond);
                        }
                    }
                }
            }
        }
    }
    
    int numCrosslinksAdd = candidateBonds.size();
    if(numCrosslinksAdd > nCrosslinks)
        numCrosslinksAdd = nCrosslinks;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    int finalNumCrosslinksAdd = 0;
    
    for(int i = 0; i < numCrosslinksAdd; i++)
    {
        int rndIndex = (int)(candidateBonds.size() * dis(threadEngine));
        Bond newBond = candidateBonds[rndIndex];
        
        candidateBonds[rndIndex] = candidateBonds.back();
        candidateBonds.pop_back();
        
        pinfo.addBond(newBond);
        finalNumCrosslinksAdd += 1;
        //bondTagsOldestToYoungest.push_back(newTag);
        
        if(candidateBonds.size() <= 0)
            break;
    }
    
    return numCrosslinksAdd;
}

//adds up to nCrosslinks, returns number of crosslinks added
int Simulation::addCrosslinksInOrder(int nCrosslinks)
{
    if(nCrosslinks <= 0)
        return 0;
    
    // search for particles between minCrosslinkBondDist and maxCrosslinkBondDist
    // should generate this list when doing excluded volume calculations and using the grid
    std::vector <pair <double, Bond>> candidateBonds;
    // add up to nCrosslinks new crosslinks from z down
    
    #pragma omp parallel for
    for(int i = 0; i < pinfo.getNumParticles(); i++)
    {
        int particleITag = pinfo.getTagAtIndex(i);
        Coordinate posI = pinfo.getPos(particleITag);
        
        std::vector <int> tempNeighbors = pinfo.getNeighbors(particleITag);
        
        for(int neighPos = 0; neighPos < tempNeighbors.size(); neighPos++)
        {
            int j = tempNeighbors[neighPos];
            
            int particleJTag = pinfo.getTagAtIndex(j);
            
            // only want to check pairs between particles once
            if(particleITag < particleJTag)
            {
                Coordinate posJ = pinfo.getPos(particleJTag);
                
                Coordinate uij = simBox.calcDisplacement(posI, posJ);
            
                double distMag = uij.getMagnitude();
                
                if(distMag >= minCrosslinkBondDist && distMag <= maxCrosslinkBondDist)
                {
                    Bond newBond = {particleITag, particleJTag, 2};
                    newBond.creationTime = currSimTime;
                    
                    // checkBondExists checks if any bond between particleITag and particleJTag exists
                    bool bondAlreadyExists = pinfo.checkBondExists(newBond);
                    
                    // no actin bond either (just define a new function checkAnyBondExists that disregards type
                    //Bond newBondTwo = {particleITag, particleJTag, 0};
                    //bool bondAlreadyExistsTwo = pinfo.checkBondExists(newBondTwo);
                    
                    if(bondAlreadyExists == false)
                    {
                        double zpos = 0.5*uij.z + posJ.z;
                        
                        #pragma omp critical
                        {
                            candidateBonds.push_back(std::make_pair(zpos, newBond));
                        }
                    }
                }
            }
        }
    }
 
    //sort candidateBonds by z value
    std::sort(candidateBonds.begin(), candidateBonds.end());
    
    int numCrosslinksAdd = candidateBonds.size();
    if(numCrosslinksAdd > nCrosslinks)
        numCrosslinksAdd = nCrosslinks;
    
    for(int i = 0; i < numCrosslinksAdd; i++)
    {
        pinfo.addBond(candidateBonds[candidateBonds.size()-1-i].second);
    }
    
    return numCrosslinksAdd;
}

/*void Simulation::calcLeadingEdgeForce(int f)
{
    int numParticles = filaments[f].getNumParticles();
        
    double totalForceMag = 0.0;
    Coordinate forceDir = Coordinate {0.0, 0.0, -1.0};
    
    double numExtraMonomers = 0.0;
    
    for(int p = 0; p < numParticles; p++)
    {
        int tempTag = filaments[f].getTagAtIndex(p);
        
        int tempPinfoPos = pinfo.getIndexOfTag(tempTag);
        
        Coordinate tempPos = pinfo.getPos(tempTag);
        
        double distToMembrane = tempPos.z;
        if(tempPos.z > -1.0)
        {
            // now caculate intersection point with p and p-1
            if(p > 0)
            {
                Coordinate disp = pinfo.getPos(filaments[f].getTagAtIndex(p-1)) - tempPos;
                disp = disp.getUnitCoord();
                
                // vertical line to boundary
                double k = tempPos.z+1.0;
                
                double length = (k*disp).getMagnitude();
                
                numExtraMonomers = length / 2.7e-3;
            }
            
            double numMonomersAboveMembrane = numExtraMonomers + (numParticles - p - 1)*14;  
            totalForceMag = numMonomersAboveMembrane * membraneForcePerMonomer;

            break;
        }
    }
    
    if(totalForceMag > 0.0)
    {
        totalForceMag = totalForceMag / (double)numParticles;
        
        for(int p = 0; p < numParticles; p++)
        {
            int tempTag = filaments[f].getTagAtIndex(p);
            pinfo.addForce(tempTag, totalForceMag*forceDir);
        }
    }
}*/

// calculates the force on each filament due to the leading edge, if a particle in a filament is above -0.5 then the entire filament
// experiences a pushing force of 10 pN distributed to each bead in the direction of the line connecting the two beads with the highest value of z
double Simulation::calcLeadingEdgeForce(int f)
{
    int numParticles = filaments[f].getNumParticles();
    
    if(numParticles < 2)
    {
        // this should not happen since this means the filament only has a single particle with no bond
        cout << "filament was smaller than size 2" << endl;
        return 0.0;
    }
        
    double totalForceMag = 0.0;
    Coordinate forceDir = Coordinate {0.0, 0.0, -1.0};
    
    double numExtraMonomers = 0.0;
    
    for(int p = 0; p < numParticles; p++)
    {
        int tempTag = filaments[f].getTagAtIndex(p);
        
        int tempPinfoPos = pinfo.getIndexOfTag(tempTag);
        
        Coordinate tempPos = pinfo.getPos(tempTag);
        
        double distToMembrane = tempPos.z;
        if(tempPos.z > -0.5)
        {
            // 1.0 pN per filament
            totalForceMag = 1.0*15;

            forceDir = pinfo.getPos(filaments[f].getTagAtIndex(numParticles-2)) - pinfo.getPos(filaments[f].getTagAtIndex(numParticles-1));
            
            if(forceDir.z > 0.0)
            {
                forceDir = -forceDir;
            }
            
            break;
        }
    }
    
    if(totalForceMag > 0.0)
    {
        Coordinate forcePerBead = totalForceMag / (double)numParticles * forceDir;
        
        for(int p = 0; p < numParticles; p++)
        {
            int tempTag = filaments[f].getTagAtIndex(p);
            pinfo.addForce(tempTag, forcePerBead);
        }
    }
    
    return totalForceMag;
}  

void Simulation::writeXYZFile(ofstream& outputXYZ)
{
    outputXYZ << pinfo.getNumParticles() << endl;
    outputXYZ << "t=" << currSimTime << endl;
    
    for(int f = 0; f < filaments.size(); f++)
    {
        int numParticlesInFilament = filaments[f].getNumParticles();
        for(int f2 = 0; f2 < numParticlesInFilament; f2++)
        {
            int p = pinfo.getIndexOfTag(filaments[f].getTagAtIndex(f2));
            Coordinate tempCoordinate = pinfo.getPosByIndex(p);
            
            Coordinate tempForce = pinfo.getForceByIndex(p); 
        
            outputXYZ << pinfo.getBeadToFilamentByIndex(p)+ 1 << " " << pinfo.getTagAtIndex(p) << " " << tempCoordinate.x << " " << tempCoordinate.y << " " << tempCoordinate.z << " " << tempForce.x << " " << tempForce.y << " " << tempForce.z;
            
            for(int nStress = 0; nStress < 2; nStress++)
            {
                // negate here so sign is the same as LAMMPS
                // 0 is bond, 1 is bending
                Coordinate tempStress = -pinfo.getStressByIndex(p,nStress);
                
                outputXYZ << " " << tempStress.x << " " << tempStress.y << " " << tempStress.z;
            }
            
            outputXYZ << " " << pinfo.getViscosityTypeByIndex(p);

            outputXYZ << endl;
        }
    }
}

void Simulation::writeBondFile(ofstream& outputBND)
{
    int numBonds = pinfo.getNumBonds();
    
    //if(numBonds > 0)
    {
        outputBND << numBonds << endl;
        outputBND << "t=" << currSimTime << endl;
        for(int i = 0; i < numBonds; i++)
        {
            Bond tempBond = pinfo.getBondByIndex(i);
            outputBND << tempBond.bondTypeIndex << " " << tempBond.i << " " << tempBond.j << " " << tempBond.creationTime << " " << tempBond.destructionTime << endl;
        }
        
        outputBND << endl;
    }
}

int Simulation::addParticlePutInGrid(Particle p, int numStressMeasurements)
{
    int newTag = pinfo.addParticle(p, numStressMeasurements);
    
    int newGrid = particleGrid.putInGrid(newTag, pinfo.getPos(newTag));
    //if(newGrid == -1)
    //    removeParticleByIndex(i);
    
    std::vector <int> tempNeighbors = particleGrid.getNeighboringTags(newTag);
    for(int i = 0; i < tempNeighbors.size(); i++)
    {
        pinfo.addNeighbor(tempNeighbors[i], newTag);
    }
    
    // also need to set neighbors of newTag as well!
    pinfo.setNeighbors(newTag, tempNeighbors);
    
    // set bond neighbors as well here???
    
    return newTag;
}

void Simulation::addFilament(int numBeads, Coordinate startPos, Coordinate direction, bool crosslinkNeighboring)
{   
    Filament newFil;
    newFil.setIsDaughter(false);
    filaments.push_back(newFil);
    int newFilamentTag = f_TaggedVector.add();
        
    int numStressMeasurements = 2;
    
    bool filamentBondedWithNeighbor = false;
    int numBondsWithNeighbors = 0;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    for(int f = 0; f < numBeads; f++)
    {
        // should put new particle in grid and update neighbors of all associated particles
        Coordinate newPosition = {startPos.x + direction.x*L_zero*f, startPos.y + direction.y*L_zero*f, startPos.z + direction.z*L_zero*f};
        newPosition = simBox.periodicWrap(newPosition);
        Particle newParticle(newPosition.x, newPosition.y, newPosition.z);
        newParticle.setType(newFilamentTag);
        
        // new bond is not incorporated into bond neighbors though!!!
        int newTag = addParticlePutInGrid(newParticle, numStressMeasurements);
        int priorTag = filaments[filaments.size()-1].addParticleBack(newTag);
        
        if(filamentBondsInBondList && priorTag > -1)
        {
            // adds a bond between priorTag and newTag of type 0 if they are on the same filament
            // as indicated from the type column in the xyz file
            Bond newBond = {priorTag, newTag, 0};
            newBond.creationTime = currSimTime;
            
            newBond.destructionTime = currSimTime + 125.0 - 5.0 * log(dis(threadEngine));
            //newBond.destructionTime = 1000000.0;
            
            //cout << "bondlifetime " << newBond.destructionTime  - currSimTime << endl;
            
            pinfo.addBond(newBond);
        }
        
        // bond to existing filament beads that are within the correct range
        // bond if appropriate to existing particles
        if(crosslinkNeighboring == true && numBondsWithNeighbors < 5)
        {
            Coordinate particleJ = pinfo.getPosByIndex(pinfo.getNumParticles()-1);
            int tempBondJTag = newTag;
            
            //std::vector <int> neighboringTags = simBox.getNeighboringTags(newTag);
            std::vector <int> neighboringTags = pinfo.getNeighbors(newTag);
        
            // this should use the grid rather than an exhaustive search
            for(int loopI = 0; loopI < neighboringTags.size(); loopI++)
            {
                int tempBondITag = neighboringTags[loopI];
                Coordinate particleI = pinfo.getPos(tempBondITag);

                Coordinate uij = particleI - particleJ;
                
                uij = simBox.periodicWrap(uij);
            
                double distMag = uij.getMagnitude();
                
                if(distMag >= minCrosslinkBondDist && distMag <= maxCrosslinkBondDist)
                {
                    //std::mt19937& threadEngine = generators[0];
                    //double randomVal = dis(threadEngine);
                    
                    // peramanent bond
                    //if(randomVal > 0.5)
                    {
                        // assumes that crosslinker bonds are stored as type 1 in bondTypes
                        Bond newBond = {tempBondITag, tempBondJTag, 1};
                        newBond.creationTime = currSimTime;
                        // never destruct permanent bonds
                        newBond.destructionTime = currSimTime + 100000.0;
                        
                        bool bondAlreadyExists = pinfo.checkBondExists(newBond);
                        
                        if(bondAlreadyExists == false)
                        {
                            pinfo.addBond(newBond);
                            numBondsWithNeighbors++;
                            //only add a single bond between new filament and existing network if applicable (arp2/3)
                            //filamentBondedWithNeighbor = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    /*cout << "F: " << newFilamentTag << endl;
    for(int f = 0; f < filaments[filaments.size()-1].getNumParticles(); f++)
    {
        int tag = filaments[filaments.size()-1].getTagAtIndex(f);
        
        std::vector <HalfBond> tempBondListB = pinfo.getBondTags(tag);
        cout << "newtag " << tag << ": ";
        for(int s = 0; s < tempBondListB.size(); s++)
        {
            cout << tempBondListB[s].j << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

void Simulation::addRandomFilament(bool crosslinkNeighboring)
{
    double randX, randY, randZ;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    randY = dis(threadEngine)*simBox.getBoxLength(1) - 0.5*simBox.getBoxLength(1);
    
    double theta = pi*0.5;
    
    randX = dis(threadEngine)*simBox.getBoxLength(0) - 0.5*simBox.getBoxLength(0);
    randZ = -1.0;
    
    double randomAngle = theta + 2.44346*dis(threadEngine) - 1.22173;
    double randomOutOfPlaneAngle = 10.0*pi/180.0 * 2.0 * (dis(threadEngine) - 0.5);
        
    double dirX = cos(randomOutOfPlaneAngle)*cos(randomAngle);
    double dirZ = cos(randomOutOfPlaneAngle)*sin(randomAngle);
    double dirY = sin(randomOutOfPlaneAngle);
    
    if(dirZ < 0.0)
    {
        dirX = copysign(1.0, dirX);
        dirZ = 0.0;
    }
    
    int filamentLength = (int)round(1.0 / L_zero) + 1;
    
    Coordinate firstBead = {randX, randY, randZ};
    Coordinate direction = {dirX, dirY, dirZ};
    
    addFilament(filamentLength, firstBead, direction, crosslinkNeighboring);
}

void Simulation::removeFilament(int filamentTag)
{
    int filamentIndex = f_TaggedVector.getIndexOfTag(filamentTag);
    
    int initialNumParticles = filaments[filamentIndex].getNumParticles();
    for(int i = 0; i < initialNumParticles; i++)
    {        
        int particleTag = filaments[filamentIndex].getBarbedTag();
        filaments[filamentIndex].removeParticleBack();
        
        pinfo.removeParticle(particleTag);
        particleGrid.removeFromGrid(particleTag);
    }
    
    filaments[filamentIndex] = filaments.back();
    filaments.pop_back();
    f_TaggedVector.remove(filamentTag);
}

void Simulation::removeFilamentsNearPos(Coordinate pos, double dist)
{
    int initialNumParticles = pinfo.getNumParticles();
    std::vector <int> removeFilamentsVector;
    
    // if a particle on a filament is near the focal adhesion mark it for removal via removeFilament function
    // should only search through particles which are near the focal adhesion
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        Coordinate currPos = pinfo.getPosByIndex(p);
        
        if(currPos.z <= pos.z + dist && currPos.z >= pos.z - dist && currPos.y <= 0.0)
        {
            Coordinate currDisp = simBox.calcDisplacement(currPos, pos);
            currPos.y = 0.0;
            
            double currDist = currDisp.getMagnitude();
            if(currDist <= dist)
            {
                int currFilamentTag = pinfo.getBeadToFilamentByIndex(p);
                
                #pragma omp critical
                {
                    bool found = false;
                    
                    for(int m = 0; m < removeFilamentsVector.size(); m++)
                    {
                        if(removeFilamentsVector[m] == currFilamentTag)
                            found = true;
                    }
                    
                    if(found == false)
                        removeFilamentsVector.push_back(currFilamentTag);
                }
            }
        }
    }
    
    for(int f = 0; f < removeFilamentsVector.size(); f++)
    {
        removeFilament(removeFilamentsVector[f]);
    }
}

void Simulation::applyNascentFAForceToTags(std::vector <int> tagList)
{
    double nascentAdhesionForceFactor = 0.1;
    
    //#pragma omp parallel for
    /*for(int p = 0; p < tagList.size(); p++)
    {
        int currTag = tagList[p];
        
        Coordinate currPos = pinfo.getPos(currTag);
        
        double realY = currPos.y;
        
        // nascent focal adhesion is only acting on particles in the bottom half of the box in the y-direction
        //if(realY < 0.0)
        {
            currPos.y = 0.0;
            Coordinate currDisp = simBox.calcDisplacement(currPos, pos);
            
            double currDist = currDisp.getMagnitude();
            if(currDist <= dist)
            {
                pinfo.setForceByIndex(p, nascentAdhesionForceFactor*pinfo.getForceByIndex(p));
            }
        }
    }*/
}

void Simulation::applyNascentFAForce(Coordinate pos, double dist)
{
   /*
    // should be able to call this only once if FA position and size are fixed
    std::vector <int> binsInFA = getBinsInRegion(pos, dist);
    // if this list is kept for multiple steps may need to deal with the particles no longer being in the list (depends on how fast these two functions are...)
    std::vector <int> tagList = getTagInBins(binsInFA);
    
    applyNascentFAForceToTags(tagList)
    */
        
    
    // if a particle on a filament is near the focal adhesion mark it for removal via removeFilament function
    // should only search through particles which are near the focal adhesion
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        Coordinate currPos = pinfo.getPosByIndex(p);
        
        if(currPos.z <= pos.z + dist && currPos.z >= pos.z - dist && currPos.y <= 0.0)
        {
            Coordinate currDisp = simBox.calcDisplacement(currPos, pos);
            currPos.y = 0.0;
            
            double currDist = currDisp.getMagnitude();
            if(currDist <= dist)
            {
                // in nascent fa
                pinfo.setForceByIndex(p, zeta / zetaInNascentFA *pinfo.getForceByIndex(p));
            }
        }
    }
}


void Simulation::removeAllCrosslinks()
{
    for(int b = pinfo.getNumBonds()-1; b > -1; b--)
    {
        Bond tempBond = pinfo.getBondByIndex(b);
        
        // crosslinks have a type of 1 or 2
        if(tempBond.bondTypeIndex != 0)
        {
            pinfo.removeBondByIndex(b);
        }
    }
}

void Simulation::removeTemporaryCrosslinks()
{
    for(int b = pinfo.getNumBonds()-1; b > -1; b--)
    {
        Bond tempBond = pinfo.getBondByIndex(b);
        
        // temporary crosslinks are type 2
        if(tempBond.bondTypeIndex == 2)
        {
            pinfo.removeBondByIndex(b);
        }
    }
}

// what happens if the bond removed is internal to a filament??? Would need to reorganize filaments like in removeParticleByIndex
void Simulation::removeBond(int bondTag)
{
    int bondIndex = pinfo.getBondPosOfTag(bondTag);
    this->removeBondByIndex(bondIndex);
}

void Simulation::removeBondByIndex(int bondIndex)
{
    Bond currBond = pinfo.getBondByIndex(bondIndex);
    
    int particleITag = currBond.i;
    int particleJTag = currBond.j;
    
    int filamentI = pinfo.getBeadToFilament(particleITag);
    int filamentJ = pinfo.getBeadToFilament(particleJTag);
    
    // remove bond before doing stuff with filaments?
    // dont run into problem with trying to delete bond after it may have already been removed (if filament of single particle was destroyed)
    // bondIndex would be the wrong index potentially if a bead was deleted and the associated bond was deleted
    pinfo.removeBondByIndex(bondIndex);
    
   // cout << "removeBondByIndex: " << filamentI << " " << filamentJ << endl;
    if(filamentI == filamentJ)
    {
        // this bond is internal to a filament, need to break this filament at this bond
        
        int filamentIndexI = f_TaggedVector.getIndexOfTag(filamentI);
        
        int indexI = filaments[filamentIndexI].getIndexOfTag(particleITag);
        int indexJ = filaments[filamentIndexI].getIndexOfTag(particleJTag);
        
        if(indexI > indexJ)
        {
            int tempSwap = indexI;
            indexI = indexJ;
            indexJ = tempSwap;
        }
        
       // cout << "break filament at " << indexI << endl;
       // cout << "break between: " << particleITag << " " << particleJTag << endl;
        
        Filament tempFilament = filaments[filamentIndexI].breakFilamentAtPos(indexI);
        
        cleanExistingFilament(filamentIndexI, filamentI);
        cleanNewFilament(tempFilament);
    }
}

// can probably do this when the grid is constructed instead (if particle is not in the grid then remove it)
void Simulation::removeParticlesOutsideBox()
{
    double minZ = simBox.getBoxLength(2);

    for(int p = pinfo.getNumParticles()-1; p > -1; p--)
    {
        int tempTag = pinfo.getTagAtIndex(p);
        
        if(tempTag >= 0 && pinfo.getPosByIndex(p).z < -minZ)
        {
            // can remove two particles at once, how to deal with this???
            // may be okay just rechecks particles which were already checked
            removeParticleByIndex(p);
        }
    }
}

void Simulation::cleanExistingFilament(int filamentIndex, int filamentTag)
{
    if(filaments[filamentIndex].getNumParticles() == 1)
    {
        // only one particle in filament, remove it and the filament
        //cout << "cleanExisting 1" << endl;
        int onlyTag = filaments[filamentIndex].getTagAtIndex(0);
        
        // removes any bonds associated with this particle as well in ParticleInfo.cpp
        pinfo.removeParticle(onlyTag);
        particleGrid.removeFromGrid(onlyTag);

        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
    }
    else if(filaments[filamentIndex].getNumParticles() == 0)
    {
        // no particles in filament, remove filament
        //cout << "cleanExisting 0" << endl;
        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
    }
}

void Simulation::cleanNewFilament(Filament tempFil)
{
    if(tempFil.getNumParticles() >= 2)
    {
        //cout << "cleanNew 2" << endl;
        
        filaments.push_back(tempFil);
        int newFilamentTag = f_TaggedVector.add();
        
        // update which filament beads are in for the new filament
        for(int v = 0; v < tempFil.getNumParticles(); v++)
        {
            int currTag = tempFil.getTagAtIndex(v);
            
            pinfo.setBeadToFilament(currTag, newFilamentTag);
        }
    }
    else if(tempFil.getNumParticles() == 1)
    {
        //cout << "cleanNew 1" << endl;
        
        int onlyTag = tempFil.getTagAtIndex(0);
        
        // always remove regardless of any bonds to onlyTag
        pinfo.removeParticle(onlyTag);
        particleGrid.removeFromGrid(onlyTag);
    }
}


void Simulation::removeParticleByIndex(int pIndex)
{
    int pTag = pinfo.getTagAtIndex(pIndex);
    
    // need to update filament indexing
    int filamentTag = pinfo.getBeadToFilamentByIndex(pIndex);
    
    Coordinate tempPos = pinfo.getPos(pTag);

    int filamentIndex = f_TaggedVector.getIndexOfTag(filamentTag);
    
    int filamentSize = filaments[filamentIndex].getNumParticles();
    
    if(filamentSize > 2)
    {        
        Filament tempFil = filaments[filamentIndex].removeBeadByTag(pTag);
        cleanExistingFilament(filamentIndex, filamentTag);
        // adds tempFil to filaments if it is big enough
        cleanNewFilament(tempFil);
    }
    else if(filamentSize == 2)
    {
        // filament with bead to be removed it only size 2, remove both bead and filament
        
        //cout << "Case 6" << endl;
        
        int otherTag = filaments[filamentIndex].getTagAtIndex(0);
        
        if(otherTag == pTag)
        {
            otherTag = filaments[filamentIndex].getTagAtIndex(1);
        }
        
        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
        
        // if this tag is of a lower index it may have not been given a bin in the grid yet by putAllParticlesInGrid
        pinfo.removeParticle(otherTag);
        particleGrid.removeFromGrid(otherTag);
    }
    else
    {
        // filament with bead to be removed is size 1 (or 0 if there is an error), remove filament
        
        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
    }
    
    
    // actually remove the particle from pinfo (also removes the bonds (but not yet angles))
    pinfo.removeParticle(pTag);
    particleGrid.removeFromGrid(pTag);
}

void Simulation::run()
{
    bool verbose = false;
    
    ofstream outputXYZ;
    outputXYZ.open ("out-" + outputFileName + ".xyz");
    
    ofstream outputBND;
    outputBND.open ("out-" + outputFileName + ".bnd");
    
    ofstream leadingEdgeFile;
    leadingEdgeFile.open ("leadingEdgeForce-" + outputFileName + ".txt");
    
    clock_t t1,t2;
	t1=clock();
    
    double initialSimTime = currSimTime;
    
    //std::cout << std::setprecision(3);
    
    omp_set_num_threads(numThreads);

    for(int i = 0; i < omp_get_max_threads(); i++)
    {
        vector <Coordinate> tempVec;
        
        for(int j = 0; j < pinfo.getNumParticles(); j++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            tempVec.push_back(c);
        }
        tempForces.push_back(tempVec);
        
        tempStresses.push_back(tempVec);
    }
    
    double desiredRetrogradeVelocity = 30.0 / 1000.0;
    double delVolume = desiredRetrogradeVelocity * simBox.getBoxLength(0) * simBox.getBoxLength(1);
    
    //double timeBetweenAddFilaments = 0.005769;
    double timeBetweenAddFilaments = numMonomersPerBead * 10 / 4.8176e5 / delVolume;
    int timestepAddFilaments = (int)round(timeBetweenAddFilaments / dt);
    cout << "timestepAddFilaments: " << timestepAddFilaments << endl;
    
    double addTemporaryCrosslinksTime = 0.025;
    int addTemporaryCrosslinksStep = (int)round(addTemporaryCrosslinksTime / dt);
    
    // temporary crosslinker density in uM
    double desiredCrosslinkingDensity = 1.0;
    int desiredNumTemporaryCrosslinks = (int)round(desiredCrosslinkingDensity * 6.022e23 * 1e-6 / 1e15 * simBox.getBoxLength(0) * simBox.getBoxLength(1) * simBox.getBoxLength(2));
    cout << "Desired num crosslinks: " << desiredNumTemporaryCrosslinks << endl;
    
    // in seconds
    double tempCrosslinkLifetime = 1.0;
    
    for(int mainLoopVar = 0; mainLoopVar < numSteps; mainLoopVar++)
    {
        if(verbose == true)
            cout << "before putting in grid" << endl;
        
        if(mainLoopVar % updateGridStep == 0)
        {
            // remove particles here if they are outside the simulation box when check if they are in grid?
            // also generate the vector for use with addCrosslinks function
            
            putAllParticlesInGrid();
        }
        
        // calculate main forces
        double leadingEdgeForce = calcForces();

        // adjust forces which are in the nascent fa
        //if(zetaInNascentFA != zeta)
        //    applyNascentFAForce(FABead, FABeadRadius);
        
        //https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
        //auto start = std::chrono::high_resolution_clock::now();
        
        moveParticles();
        
        currSimTime += dt;
        
        if(mainLoopVar % 100 == 0)
            leadingEdgeFile << currSimTime << " " << leadingEdgeForce << endl;
        
        
        // update actin network connectivity here
        
        
        int numRemovedCrosslinks = removeFlaggedBonds();
        

        // add a random filament to the network if it is time to do so
        if(mainLoopVar % timestepAddFilaments == 0)
        {
            addRandomFilament(crosslinkNeighboringUponAdding);
        }
        
        // remove filaments that have a bead within FABeadRadius from FABead
        // removeFilamentsNearPos(FABead, FABeadRadius);
        
        
        // add crosslinks up to desiredNumCrosslinks
        if(mainLoopVar % addTemporaryCrosslinksStep == 0)
        {
            int currNumTempCrosslinks = 0;
            
            // need to count the number of temporary crosslinks here since the crosslinks could have been removed from other means than just aging
            // would need to have a method of accounting for removal of all bonds, this may be possible through removeBondByIndex wrapper function in this class
            for(int b = 0; b < pinfo.getNumBonds(); b++)
            {
                if(pinfo.getBondByIndex(b).bondTypeIndex == 2)
                {
                    currNumTempCrosslinks++;
                }
            }
            
            if(currNumTempCrosslinks < desiredNumTemporaryCrosslinks)
            {
                int numAddedCrosslinks = addTemporaryCrosslinks(desiredNumTemporaryCrosslinks - currNumTempCrosslinks);
                
                if(verbose)
                    cout << "Added crosslinkers: " << numAddedCrosslinks << ", total num crosslinkers: " << numAddedCrosslinks + currNumTempCrosslinks << endl;
            }
        }
        
        if(verbose == true)
            cout << "before snapshot" << endl;
        
        // take a snapshot
        if(mainLoopVar % snapshotStep == 0)
        {
            // remove particles outside of the box for the snapshot            
            removeParticlesOutsideBox();
            
            cout << fixed << currSimTime << " / " << (totalSimulationTime + initialSimTime) << endl;
            writeXYZFile(outputXYZ);
            writeBondFile(outputBND);
        }       
    }
    
    t2=clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;

	cout<<seconds<<endl;
    
    outputXYZ.close();
    outputBND.close();
}