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

#ifndef TAGGEDVECTOR_H
#define TAGGEDVECTOR_H

#include <vector>
#include <queue>
#include <stdexcept>

class TaggedVector
{
    private:
        unsigned int numCurrEntries;
        unsigned int maxHistoricalTag;
        std::vector <int> v_tag;
        std::vector <int> v_rtag;
        std::queue <int> q_recycledTags;
        
    public:
        TaggedVector();
        
        int getNumEntries();
        
        int add();
		void remove(int oldTag);
        
        int getTagAtIndex(int pos);
		int getIndexOfTag(int tag);
        
        int getCurrMaxTag();
        
        void setTagAtPos(int pos, int newTag);
};

#endif