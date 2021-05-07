#ifndef QUEST_CACHE_H
#define QUEST_CACHE_H
#include "QuEST.h"
#include "QuEST_gates.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

const int CACHE_WIDTH = 16;
const int BATCH_SIZE = 8;
const int64_t CACHE_SIZE = 1LL << CACHE_WIDTH;
const int64_t numTasks = CACHE_SIZE >> 1;

#define defineHalfBlockSize(g, halfBlock) \
    const int64_t halfBlock = (1LL << g->targetQubit) - 1

#define defineIndices(halfBlock, thisTask, indexUp, indexLo) \
    int64_t indexUp = (thisTask & ~halfBlock)*2LL + (thisTask & halfBlock); \
    int64_t indexLo = indexUp + halfBlock + 1

#define getControlBit(idx, g) \
    ((idx + g->chunkOffset) >> g->controlQubit) & 1

#define controlBlockSize(g) (1LL << g->controlQubit)

//------------------------------------ copy some tool functions
static int chunkIsUpper(int chunkId, long long int chunkSize, int targetQubit)
{       
    long long int sizeHalfBlock = 1LL << (targetQubit);
    long long int sizeBlock = sizeHalfBlock*2;
    long long int posInBlock = (chunkId*chunkSize) % sizeBlock;
    return posInBlock<sizeHalfBlock;
}
static int getChunkPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit)
{
    long long int sizeHalfBlock = 1LL << (targetQubit);
    int chunksPerHalfBlock = sizeHalfBlock/chunkSize;
    if (chunkIsUpper){
        return chunkId + chunksPerHalfBlock;
    } else {
        return chunkId - chunksPerHalfBlock;
    }
}
static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}
//------------------------------------ copy some tool functions


static void hadamardCache(GateObject *g, const int64_t cacheid)
{
    defineHalfBlockSize(g, halfBlock);
    qreal recRoot2 = 1.0/sqrt(2);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];

        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];

        g->svreal[indexUp] = recRoot2*(stateRealUp + stateRealLo);
        g->svimag[indexUp] = recRoot2*(stateImagUp + stateImagLo);

        g->svreal[indexLo] = recRoot2*(stateRealUp - stateRealLo);
        g->svimag[indexLo] = recRoot2*(stateImagUp - stateImagLo);
    }
}


static void tGateCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
    const qreal angle = 1/sqrt(2);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];
        g->svreal[indexLo] = angle*stateRealLo - angle*stateImagLo;
        g->svimag[indexLo] = angle*stateRealLo + angle*stateImagLo;
    }
}


static void sGateCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
	volatile qreal cosAngle = 0;
	volatile qreal sinAngle = 1;
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];
        g->svreal[indexLo] = cosAngle*stateRealLo - sinAngle*stateImagLo;
        g->svimag[indexLo] = sinAngle*stateRealLo + cosAngle*stateImagLo;
    }
}


static void pauliXCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];

        g->svreal[indexUp] = g->svreal[indexLo];
        g->svimag[indexUp] = g->svimag[indexLo];

        g->svreal[indexLo] = stateRealUp;
        g->svimag[indexLo] = stateImagUp;
    }
}


static void pauliYCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];

        g->svreal[indexUp] = g->svimag[indexLo];
        g->svimag[indexUp] = -g->svreal[indexLo];
        g->svreal[indexLo] = -stateImagUp;
        g->svimag[indexLo] = stateRealUp;
    }
}


static void pauliZCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
	volatile qreal cosAngle = -1;
	volatile qreal sinAngle = 0;
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];
        g->svreal[indexLo] = cosAngle*stateRealLo - sinAngle*stateImagLo;
        g->svimag[indexLo] = sinAngle*stateRealLo + cosAngle*stateImagLo;
    }
}


static void rotateXCache(GateObject *g, const int64_t cacheid) {
	Vector unitAxis = {1, 0, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;
    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        // store current state vector values in temp variables
        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];
        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];

        // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
        g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
        g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

        // state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
        g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
        g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
    }
}


static void rotateYCache(GateObject *g, const int64_t cacheid) {
	Vector unitAxis = {0, 1, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        // store current state vector values in temp variables
        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];
        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];

        // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
        g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
        g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

        // state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
        g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
        g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
    }
}


static void rotateZCache(GateObject *g, const int64_t cacheid) {
	Vector unitAxis = {0, 0, 1};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        // store current state vector values in temp variables
        qreal stateRealUp = g->svreal[indexUp];
        qreal stateImagUp = g->svimag[indexUp];

        qreal stateRealLo = g->svreal[indexLo];
        qreal stateImagLo = g->svimag[indexLo];

        // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
        g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
        g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

        // state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
        g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
        g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
    }
}


static void controlledNotCache(GateObject *g, const int64_t cacheid)
{
    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        int controlBit = getControlBit(indexUp, g);
        if (controlBit){
            qreal stateRealUp = g->svreal[indexUp];
            qreal stateImagUp = g->svimag[indexUp];

            g->svreal[indexUp] = g->svreal[indexLo];
            g->svimag[indexUp] = g->svimag[indexLo];

            g->svreal[indexLo] = stateRealUp;
            g->svimag[indexLo] = stateImagUp;
        }
    }
}


static void controlledPauliYCache(GateObject *g, const int64_t cacheid) {
    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);

        if (getControlBit(indexUp, g)) {
            qreal stateRealUp = g->svreal[indexUp];
            qreal stateImagUp = g->svimag[indexUp];

            // update under +-{{0, -i}, {i, 0}}
            g->svreal[indexUp] = g->svimag[indexLo];
            g->svimag[indexUp] = -g->svreal[indexLo];
            g->svreal[indexLo] = -stateImagUp;
            g->svimag[indexLo] = stateRealUp;
        }
    }
}


static void controlledRotateXCache(GateObject *g, const int64_t cacheid) {
	Vector unitAxis = {1, 0, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);
        
        if (getControlBit(indexUp, g)) {
            // store current state vector values in temp variables
            qreal stateRealUp = g->svreal[indexUp];
            qreal stateImagUp = g->svimag[indexUp];
            qreal stateRealLo = g->svreal[indexLo];
            qreal stateImagLo = g->svimag[indexLo];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
			g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
			g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
        }
    }
}


static void controlledRotateYCache(GateObject *g, const int64_t cacheid) {
	Vector unitAxis = {0, 1, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

    defineHalfBlockSize(g, halfBlock);
    for (int64_t t0 = 0; t0 < numTasks; t0++) {
        int64_t thisTask = t0 + cacheid * numTasks;
        defineIndices(halfBlock, thisTask, indexUp, indexLo);
        
        if (getControlBit(indexUp, g)) {
            // store current state vector values in temp variables
            qreal stateRealUp = g->svreal[indexUp];
            qreal stateImagUp = g->svimag[indexUp];
            qreal stateRealLo = g->svreal[indexLo];
            qreal stateImagLo = g->svimag[indexLo];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
			g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
			g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
        }
    }
}


static void controlledRotateZCache(GateObject *g, const int64_t cacheid)
{
	Vector unitAxis = {0, 0, 1};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;
    defineHalfBlockSize(g, halfBlock);
#define KERNEL(g, indexUp, indexLo) \
    qreal stateRealUp = g->svreal[indexUp]; \
    qreal stateImagUp = g->svimag[indexUp]; \
    qreal stateRealLo = g->svreal[indexLo]; \
    qreal stateImagLo = g->svimag[indexLo]; \
	g->svreal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo; \
	g->svimag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo; \
	g->svreal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo; \
	g->svimag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
    if (!BATCH_SIZE || halfBlock < BATCH_SIZE || controlBlockSize(g) < BATCH_SIZE) {
        for (int64_t t0 = 0; t0 < numTasks; t0++) {
            int64_t thisTask = t0 + cacheid * numTasks;
            defineIndices(halfBlock, thisTask, indexUp, indexLo);
            if (getControlBit(indexUp, g)) {
                KERNEL(g, indexUp, indexLo);
            }
        }
    } else {
        for (int64_t t0 = 0; t0 < numTasks; t0 += BATCH_SIZE) {
            int64_t thisTask = t0 + cacheid * numTasks;
            defineIndices(halfBlock, thisTask, indexUp, indexLo);
            if (getControlBit(indexUp, g)) {
                for (int i = 0; i < BATCH_SIZE; i++, indexUp++, indexLo++) {
                    KERNEL(g, indexUp, indexLo);
                }
            }
        }
    }
#undef KERNEL
}


#ifdef __cplusplus
}
#endif

#endif //QUEST_DISTRIBUTE_CACHE_H
