/*
Copyright 2022 Matthias Müller - Ten Minute Physics
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef SPATIALHASH_H
#define SPATIALHASH_H


#include <cstdlib>
#include <cmath>
#include <QVector>
#include <algorithm>
#include <iostream>

#include "particle.h"

class SpatialHash {

public:

    // Constructor - Store the input values to member variables
    SpatialHash(double spacing, unsigned int maxNumObjects)
    {
        this->maxNumObjects = maxNumObjects;
        this->spacing = spacing;
        this->tableSize = 5 * maxNumObjects;
        this->cellStart.resize(tableSize + 1);
        this->cellEntries.resize(maxNumObjects);
        this->queryIds.resize(maxNumObjects);
        this->querySize = 0;
    }

    // hash function from the slides
    unsigned int hashCoords(int xi, int yi, int zi){
        int h = (xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481); // fantasy function
        return std::abs(h) % this->tableSize;
    }

    // compute the coords of the cell that contains the object with the given coordinate
    int intCoord(double coord){
        return std::floor(coord / this->spacing);
    }

    // takes the pos of the object and returns the index of the containing cell in the hash table
    unsigned int hashPos(QVector<Particle *> particles, unsigned int nr){
        return hashCoords(
            this->intCoord(particles[nr]->pos.x()),
            this->intCoord(particles[nr]->pos.y()),
            this->intCoord(particles[nr]->pos.z())
        );
    }

    // creates the hash, given the position of all the objects
    void create(QVector<Particle *> particles){
        unsigned int numObjects = std::min(particles.size(), this->cellEntries.size());

        // determine cell sizes

        this->cellStart.fill(0);
        this->cellEntries.fill(0);

        for(unsigned int i=0; i < numObjects; i++){
            unsigned int h = this->hashPos(particles,i);
            this->cellStart[h]++;
        }

        // determine cells starts

        unsigned int start = 0;
        for(unsigned int i=0; i < this->tableSize; i++){
            start += this->cellStart[i];
            this->cellStart[i] = start;
        }
        this->cellStart[this->tableSize] = start; //guard

        // fill in objects ids

        for(unsigned int i=0; i < numObjects; i++){
            unsigned int h = this->hashPos(particles,i);
            this->cellStart[h]--;
            this->cellEntries[this->cellStart[h]] = i;
        }
    }

    // how to retrieve objects from the hash
    void query(QVector<Particle *> particles, unsigned int nr, float maxDist){
        int x0 = this->intCoord(particles[nr]->pos.x() - maxDist);
        int y0 = this->intCoord(particles[nr]->pos.y() - maxDist);
        int z0 = this->intCoord(particles[nr]->pos.z() - maxDist);

        int x1 = this->intCoord(particles[nr]->pos.x() + maxDist);
        int y1 = this->intCoord(particles[nr]->pos.y() + maxDist);
        int z1 = this->intCoord(particles[nr]->pos.z() + maxDist);

        this->querySize = 0;

        for(int  xi=x0; xi<=x1; xi++){
            for(int  yi=y0;yi<=y1; yi++){
                for(int  zi=z0;zi<=z1; zi++){
                    unsigned int h = this->hashCoords(xi,yi,zi);
                    unsigned int start = this->cellStart[h];
                    unsigned int end = this->cellStart[h+1];

                    for(unsigned int i=start; i<end; i++){
                        this->queryIds[this->querySize] = this->cellEntries[i];
                        this->querySize++;
                    }
                }
            }
        }
    }


    // ---- variables ---- //
    double spacing;
    unsigned int tableSize, querySize;
    QVector<unsigned int> cellStart, cellEntries, queryIds;
    unsigned int maxNumObjects;

};

#endif // SPATIALHASH_H
