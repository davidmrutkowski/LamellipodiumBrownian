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
 
#include "Grid.h"

using namespace std;

Grid::Grid()
{
    gridSize = Coordinate {1.0, 1.0, 1.0};
    periodicX = true;
    periodicY = true;
    periodicZ = true;
}

Grid::Grid(Coordinate b, bool px, bool py, bool pz)
{
    gridSize = b;
    
    periodicX = px;
    periodicY = py;
    periodicZ = pz;
}

void Grid::setGridDimensions(Coordinate newGridSize)
{
	gridSize = newGridSize;
}

int Grid::getNumBins()
{
	return grid.size();
}

void Grid::setNumBins(int newNumBinX, int newNumBinY, int newNumBinZ)
{
	numBinX = newNumBinX;
	numBinY = newNumBinY;
	numBinZ = newNumBinZ;
	
	grid.resize(numBinX*numBinY*numBinZ);
	
    
    
	binSize.x = gridSize.x / (double)numBinX;
	binSize.y = gridSize.y / (double)numBinY;
	binSize.z = gridSize.z / (double)numBinZ;
}

// sets binSize to the closest value that is larger possible based on desiredBinSize in each dimension
void Grid::setBinSize(double desiredBinSize)
{
    int tempNumBinX = (int)floor(gridSize.x / desiredBinSize);
    int tempNumBinY = (int)floor(gridSize.y / desiredBinSize);
    int tempNumBinZ = (int)floor(gridSize.z / desiredBinSize);
    
    
    setNumBins(tempNumBinX, tempNumBinY, tempNumBinZ);
}

void Grid::setPeriodic(bool newPeriodicX, bool newPeriodicY, bool newPeriodicZ)
{
	periodicX = newPeriodicX;
	periodicY = newPeriodicY;
	periodicZ = newPeriodicZ;
}

int Grid::putInGrid(int tag, Coordinate c)
{
	int binx = (int)((c.x + gridSize.x*0.5) / (double)binSize.x);
	int biny = (int)((c.y + gridSize.y*0.5) / (double)binSize.y);
	// z values are always negative, if not have to change this	
	int binz = (int)(-c.z / (double)binSize.z);
	
	// fix for if bead is at boundary
	if(binx == numBinX && periodicX)
	{
		binx = 0;
	}
	if(biny == numBinY && periodicY)
	{
		biny = 0;
	}
	if(binz == numBinZ && periodicZ)
	{
		binz = 0;
	}
	
	int bin = (binx*numBinY*numBinZ) + biny*numBinZ + binz;
	
	if(bin < getNumBins() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
	{
        // only resize tags if it needs to be larger to accommodate tag
        if(tag >= tags.size())
        {
            tags.resize(tag+1);
        }
        else if(tags[tag] != -1)
        {
            // this tag is already present in the grid somewhere, remove it from where it is and then add this new position
            removeFromGrid(tag);
        }

        grid[bin].push_back(tag);
        tags[tag] = bin;
        
		return bin;
	}
	else
	{
        // coordinate c is outside of grid so do not give it a position
        if(tag >= tags.size())
            tags.resize(tag+1);
        tags[tag] = -1;
        
		return -1;
	}
}

bool Grid::removeFromGrid(int tag)
{
    if(tag >= 0 && tag < tags.size())
    {
        int bin = tags[tag];
        
        // particle is not in grid but may be present in simulation
        if(bin == -1)
            return false;
        
        int binSize = grid[bin].size();
        for(int i = 0; i < binSize; i++)
        {
            if(tag == grid[bin][i])
            {
                grid[bin][i] = grid[bin].back();
                grid[bin].pop_back();
                
                tags[tag] = -1;
                
                return true;
            }
        }
        
        /*cout << tag << " " << tags[tag] << endl;
        for(int i = 0; i < grid.size(); i++)
        {
            for(int j = 0; j < grid[i].size(); j++)
            {
                if(grid[i][j] == tag)
                {
                    cout << tags.size() << endl;
                    cout << "found: " << tag << " " << i << " " << j << endl;
                    cout << "instead of in bin: " << bin << " " << tags[tag] << endl;
                }
                
            }
        }*/
    }
    
    //cout << "attempted to remove from grid and failed: " << tag << endl;
    return false;
}

void Grid::clearGrid()
{
	for(int i = 0; i < grid.size(); i++)
	{
		grid[i].clear();
	}
    
    tags.clear();
}

int Grid::getBin(Coordinate c)
{
    int binx = (int)((c.x + gridSize.x*0.5) / (double)binSize.x);
	int biny = (int)((c.y + gridSize.y*0.5) / (double)binSize.y);
	// z values are always negative, if not have to change this	
	int binz = (int)(-c.z / (double)binSize.z);
	
	// fix for if bead is at boundary
	if(binx == numBinX && periodicX)
	{
		binx = 0;
	}
	if(biny == numBinY && periodicY)
	{
		biny = 0;
	}
	if(binz == numBinZ && periodicZ)
	{
		binz = 0;
	}
	
	int bin = (binx*numBinY*numBinZ) + biny*numBinZ + binz;
	
	if(bin < getNumBins() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
    {
        return bin;
    }
    
    return -1;
}

std::tuple<int, int, int> Grid::binToGridPositions(int binIndex)
{
    int gridX = (int)(binIndex / (double)numBinY / (double)numBinZ);
	binIndex = binIndex - gridX * numBinY * numBinZ;
	int gridY = (int)(binIndex / (double)numBinZ);
	binIndex = binIndex - gridY * numBinZ;
	int gridZ = binIndex;
    
    return std::make_tuple(gridX, gridY, gridZ);
}

int Grid::gridPositionsToBin(int binx, int biny, int binz)
{
    if(periodicX)
        binx = binx - numBinX*(int)floor((double)(binx)/(double)(numBinX));
    if(periodicY)
        biny = biny - numBinY*(int)floor((double)(biny)/(double)(numBinY));
    if(periodicZ)
        binz = binz - numBinZ*(int)floor((double)(binz)/(double)(numBinZ));
    
    return binx*numBinY*numBinZ + biny*numBinZ + binz;
}

std::vector<int> Grid::getBinsInRegion(Coordinate centerPosition, double radius)
{
    std::vector <int> binList;
    
    int centralBin = getBin(centerPosition);
    
    int binx = (int)((centerPosition.x + gridSize.x*0.5) / (double)binSize.x);
	int biny = (int)((centerPosition.y + gridSize.y*0.5) / (double)binSize.y);
	// z values are always negative, if not have to change this	
	int binz = (int)(-centerPosition.z / (double)binSize.z);
	
	// fix for if bead is at boundary
	if(binx == numBinX && periodicX)
	{
		binx = 0;
	}
	if(biny == numBinY && periodicY)
	{
		biny = 0;
	}
	if(binz == numBinZ && periodicZ)
	{
		binz = 0;
	}
	
	//int centralBin = (binx*numBinY*numBinZ) + biny*numBinZ + binz;
	
    if(centralBin < getNumBins() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
    {    
        int numX = (int)ceil(radius / binSize.x);
        int numY = (int)ceil(radius / binSize.y);
        int numZ = (int)ceil(radius / binSize.z);
        
        for(int i = -numX; i <= numX; i++)
        {
            int currBinX = i + binx;
            
            for(int j = -numY; j <= numY; j++)
            {
                int currBinY = j + biny;
                
                for(int k = -numZ; k <= numZ; k++)
                {
                    int currBinZ = k + binz;
                    
                    int bin = (currBinX*numBinY*numBinZ) + currBinY*numBinZ + currBinZ;
                    
                    if(bin < getNumBins() && currBinX >= 0 && currBinX < numBinX && currBinY >= 0 && currBinY < numBinY && currBinZ >= 0 && currBinZ < numBinZ)
                    {
                        binList.push_back(bin);
                    }
                }
            }
        }
    }
    
    return binList;
}

std::vector<int> Grid::getTagInBins(std::vector <int> binList)
{
    std::vector <int> tagList;
    
    for(int b = 0; b < binList.size(); b++)
    {
        int bin = binList[b];
        for(int g = 0; g < grid[bin].size(); g++)
        {
            int newTag = grid[bin][g];
            
            tagList.push_back(newTag);                                
        }
    }
        
    return tagList;
}

std::vector<int> Grid::getTagsInRegion(Coordinate centerPosition, double radius)
{
    std::vector <int> tagList;
    
    int centralBin = getBin(centerPosition);
    
    int binx = (int)((centerPosition.x + gridSize.x*0.5) / (double)binSize.x);
	int biny = (int)((centerPosition.y + gridSize.y*0.5) / (double)binSize.y);
	// z values are always negative, if not have to change this	
	int binz = (int)(-centerPosition.z / (double)binSize.z);
	
	// fix for if bead is at boundary
	if(binx == numBinX && periodicX)
	{
		binx = 0;
	}
	if(biny == numBinY && periodicY)
	{
		biny = 0;
	}
	if(binz == numBinZ && periodicZ)
	{
		binz = 0;
	}
	
	//int centralBin = (binx*numBinY*numBinZ) + biny*numBinZ + binz;
	
    if(centralBin < getNumBins() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
    {    
        int numX = (int)ceil(radius / binSize.x);
        int numY = (int)ceil(radius / binSize.y);
        int numZ = (int)ceil(radius / binSize.z);
        
        for(int i = -numX; i <= numX; i++)
        {
            int currBinX = i + binx;
            
            for(int j = -numY; j <= numY; j++)
            {
                int currBinY = j + biny;
                
                for(int k = -numZ; k <= numZ; k++)
                {
                    int currBinZ = k + binz;
                    
                    int bin = (currBinX*numBinY*numBinZ) + currBinY*numBinZ + currBinZ;
                    
                    if(bin < getNumBins() && currBinX >= 0 && currBinX < numBinX && currBinY >= 0 && currBinY < numBinY && currBinZ >= 0 && currBinZ < numBinZ)
                    {
                        for(int g = 0; g < grid[bin].size(); g++)
                        {
                            int newTag = grid[bin][g];
                            
                            tagList.push_back(newTag);                                
                        }
                    }
                }
            }
        }
    }
    
	return tagList;
}

std::vector<int> Grid::getNeighboringTags(int particleTag)
{
	std::vector <int> neighbors;
	
    int binIndex = tags[particleTag];
    
    if(binIndex < 0)
        return neighbors;
    
	int gridX, gridY, gridZ;
    
    std::tie (gridX, gridY, gridZ) = binToGridPositions(binIndex);
	
	for(int m = -1; m <= 1; m++)
	{
		for(int n = -1; n <= 1; n++)
		{
			for(int o = -1; o <= 1; o++)
			{
				int binx = gridX + m;
				int biny = gridY + n;
				int binz = gridZ + o;
				
				int tempGrid = gridPositionsToBin(binx, biny, binz);
				
				// if the neighboring bin is within the grid then add the neighbors to neighbors
				if(tempGrid < grid.size() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
				{
					for(int g = 0; g < grid[tempGrid].size(); g++)
					{
						int newTag = grid[tempGrid][g];
						
                        if(newTag != particleTag)
                            neighbors.push_back(newTag);
					}
				}
			}
		}
	}
	
	return neighbors;
}

// finds neighbors of either particleTagI or particleTagJ without duplicates
std::vector<int> Grid::getNeighboringTags(int particleTagI, int particleTagJ)
{
	std::vector <int> neighbors;
	
    if(particleTagI >= tags.size() || particleTagJ >= tags.size())
        return neighbors;
    
    int binIndexI = tags[particleTagI];
    int binIndexJ = tags[particleTagJ];
    
    if(binIndexI < 0 || binIndexJ < 0)
        return neighbors;
    
    if(binIndexI == binIndexJ)
        return getNeighboringTags(particleTagI);
    
    int gridXI, gridYI, gridZI;
    std::tie (gridXI, gridYI, gridZI) = binToGridPositions(binIndexI);
    
    int gridXJ, gridYJ, gridZJ;
    std::tie (gridXJ, gridYJ, gridZJ) = binToGridPositions(binIndexJ);
    
    std::vector <int> binList;
    for(int m = -1; m <= 1; m++)
	{
		for(int n = -1; n <= 1; n++)
		{
			for(int o = -1; o <= 1; o++)
			{
				int binx = gridXI + m;
				int biny = gridYI + n;
				int binz = gridZI + o;
                
                int tempGrid = binx*numBinY*numBinZ + biny*numBinZ + binz;
				
				// if the neighboring bin is within the grid then add the neighbors to neighbors
				if(tempGrid < grid.size() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
				{
                    binList.push_back(tempGrid);
                }
            }
        }
    }
    for(int m = -1; m <= 1; m++)
	{
		for(int n = -1; n <= 1; n++)
		{
			for(int o = -1; o <= 1; o++)
			{
				int binx = gridXJ + m;
				int biny = gridYJ + n;
				int binz = gridZJ + o;
                
                int tempGrid = binx*numBinY*numBinZ + biny*numBinZ + binz;
				
				// if the neighboring bin is within the grid then add the neighbors to neighbors
				if(tempGrid < grid.size() && binx >= 0 && binx < numBinX && biny >= 0 && biny < numBinY && binz >= 0 && binz < numBinZ)
				{     
                    // if tempGrid was not already found in binList then add it
                    if(std::find(binList.begin(), binList.end(), tempGrid) == binList.end())
                        binList.push_back(tempGrid);
                }
            }
        }
    }
    
    for(int i = 0; i < binList.size(); i++)
    {
        int tempGrid = binList[i];
        
        for(int g = 0; g < grid[tempGrid].size(); g++)
        {
            int newTag = grid[tempGrid][g];
            
            if(newTag != particleTagI && newTag != particleTagJ)
                neighbors.push_back(newTag);
        }
    }
	
	return neighbors;
}

std::vector<int> Grid::getTagsAboveZVal(double zPos)
{
    int highestZBin = (int)(-zPos / (double)binSize.z);
    
    std::vector<int> particleTags;
    
    while(highestZBin >= 0)
    {
        for(int i = 0; i < numBinX*numBinY; i++)
        {
            int currBin = highestZBin + i*numBinZ;
            
            for(int j = 0; j < grid[currBin].size(); j++)
            {
                particleTags.push_back(grid[currBin][j]);
            }
        }
        
        highestZBin--;
    }
    
    return particleTags;
}	