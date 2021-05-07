#ifndef QUEST_DISTRIBUTE_CACHE_H
#define QUEST_DISTRIBUTE_CACHE_H
#include "QuEST.h"
#include "QuEST_gates.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void hadamardHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    qreal recRoot2 = 1.0/sqrt(2);
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealUp = stateVecRealUp[enclose_offset];
		stateImagUp = stateVecImagUp[enclose_offset];

		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		stateVecRealUp[enclose_offset] = recRoot2*(stateRealUp + stateRealLo);
		stateVecImagUp[enclose_offset] = recRoot2*(stateImagUp + stateImagLo);

		stateVecRealLo[enclose_offset] = recRoot2*(stateRealUp - stateRealLo);
		stateVecImagLo[enclose_offset] = recRoot2*(stateImagUp - stateImagLo);
	} 
}

void pauliXHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealOut,stateImagOut;
    long long int thisTask;         
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealIn=qureg.pairStateVec.real, *stateVecImagIn=qureg.pairStateVec.imag;
    qreal *stateVecRealOut=qureg.stateVec.real, *stateVecImagOut=qureg.stateVec.imag;
	if (updateUpper)
		offset = 0;
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealOut = stateVecRealOut[enclose_offset];
		stateImagOut = stateVecImagOut[enclose_offset];

		stateVecRealOut[enclose_offset] = stateVecRealIn[enclose_offset];
		stateVecImagOut[enclose_offset] = stateVecImagIn[enclose_offset];

		stateVecRealIn[enclose_offset] = stateRealOut;
		stateVecImagIn[enclose_offset] = stateImagOut;
	} 
}

void pauliYHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}
	long long int enclose_offset;
	int conjFac = 1;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealUp = stateVecRealUp[enclose_offset];
		stateImagUp = stateVecImagUp[enclose_offset];

		stateVecRealUp[enclose_offset] = conjFac * stateVecImagLo[enclose_offset];
		stateVecImagUp[enclose_offset] = conjFac * -stateVecRealLo[enclose_offset];
		stateVecRealLo[enclose_offset] = conjFac * -stateImagUp;
		stateVecImagLo[enclose_offset] = conjFac * stateRealUp;
	} 
}

void rotateXHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}
	Vector unitAxis = {1, 0, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealUp = stateVecRealUp[enclose_offset];
		stateImagUp = stateVecImagUp[enclose_offset];

		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
		stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
			- betaReal*stateRealLo - betaImag*stateImagLo;
		stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
			- betaReal*stateImagLo + betaImag*stateRealLo;

		// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
		stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
			+ alphaReal*stateRealLo + alphaImag*stateImagLo;
		stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
			+ alphaReal*stateImagLo - alphaImag*stateRealLo;
	} 
}

void rotateYHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}
	Vector unitAxis = {0, 1, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealUp = stateVecRealUp[enclose_offset];
		stateImagUp = stateVecImagUp[enclose_offset];

		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
		stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
			- betaReal*stateRealLo - betaImag*stateImagLo;
		stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
			- betaReal*stateImagLo + betaImag*stateRealLo;

		// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
		stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
			+ alphaReal*stateRealLo + alphaImag*stateImagLo;
		stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
			+ alphaReal*stateImagLo - alphaImag*stateRealLo;
	} 
}

void rotateZHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}
	Vector unitAxis = {0, 0, 1};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		stateRealUp = stateVecRealUp[enclose_offset];
		stateImagUp = stateVecImagUp[enclose_offset];

		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
		stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
			- betaReal*stateRealLo - betaImag*stateImagLo;
		stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
			- betaReal*stateImagLo + betaImag*stateRealLo;

		// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
		stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
			+ alphaReal*stateRealLo + alphaImag*stateImagLo;
		stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
			+ alphaReal*stateImagLo - alphaImag*stateRealLo;
	} 
}

void controlledNotHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
	int controlQubit = g->controlQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    int controlBit;
	long long int enclose_offset;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		controlBit = extractBit(controlQubit, enclose_offset+chunkId*chunkSize);
		if (controlBit){
			stateRealUp = stateVecRealUp[enclose_offset];
			stateImagUp = stateVecImagUp[enclose_offset];

			stateVecRealUp[enclose_offset] = stateVecRealLo[enclose_offset];
			stateVecImagUp[enclose_offset] = stateVecImagLo[enclose_offset];

			stateVecRealLo[enclose_offset] = stateRealUp;
			stateVecImagLo[enclose_offset] = stateImagUp;
		}
	} 
}

void controlledPauliYHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
	int controlQubit = g->controlQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    int controlBit;
	long long int enclose_offset;
	int conjFac = 1;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		controlBit = extractBit(controlQubit, enclose_offset+chunkId*chunkSize);
		if (controlBit){
			stateRealUp = stateVecRealUp[enclose_offset];
			stateImagUp = stateVecImagUp[enclose_offset];

			// update under +-{{0, -i}, {i, 0}}
			stateVecRealUp[enclose_offset] = conjFac * stateVecImagLo[enclose_offset];
			stateVecImagUp[enclose_offset] = conjFac * -stateVecRealLo[enclose_offset];
			stateVecRealLo[enclose_offset] = conjFac * -stateImagUp;
			stateVecImagLo[enclose_offset] = conjFac * stateRealUp;
		}
	} 
}

void controlledRotateXHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
	int controlQubit = g->controlQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    int controlBit;
	long long int enclose_offset;

	Vector unitAxis = {1, 0, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		controlBit = extractBit (controlQubit, enclose_offset+chunkId*chunkSize);
		if (controlBit){
			// store current state vector values in temp variables
			stateRealUp = stateVecRealUp[enclose_offset];
			stateImagUp = stateVecImagUp[enclose_offset];

			stateRealLo = stateVecRealLo[enclose_offset];
			stateImagLo = stateVecImagLo[enclose_offset];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
				- betaReal*stateRealLo - betaImag*stateImagLo;
			stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
				- betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
				+ alphaReal*stateRealLo + alphaImag*stateImagLo;
			stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
				+ alphaReal*stateImagLo - alphaImag*stateRealLo;
		}
	} 
} 

void controlledRotateYHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
	int controlQubit = g->controlQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    int controlBit;
	long long int enclose_offset;

	Vector unitAxis = {0, 1, 0};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		controlBit = extractBit (controlQubit, enclose_offset+chunkId*chunkSize);
		if (controlBit){
			// store current state vector values in temp variables
			stateRealUp = stateVecRealUp[enclose_offset];
			stateImagUp = stateVecImagUp[enclose_offset];

			stateRealLo = stateVecRealLo[enclose_offset];
			stateImagLo = stateVecImagLo[enclose_offset];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
				- betaReal*stateRealLo - betaImag*stateImagLo;
			stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
				- betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
				+ alphaReal*stateRealLo + alphaImag*stateImagLo;
			stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
				+ alphaReal*stateImagLo - alphaImag*stateRealLo;
		}
	} 
} 

void controlledRotateZHalf(GateObject *g, const int64_t cacheid)
{
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;
	int controlQubit = g->controlQubit;
    int updateUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
    qreal *stateVecRealUp = qureg.pairStateVec.real, *stateVecImagUp = qureg.pairStateVec.imag;
    qreal *stateVecRealLo = qureg.stateVec.real, *stateVecImagLo = qureg.stateVec.imag;
	if (updateUpper) {
		offset = 0;
		stateVecRealUp = qureg.stateVec.real;
		stateVecImagUp = qureg.stateVec.imag;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    int controlBit;
	long long int enclose_offset;

	Vector unitAxis = {0, 0, 1};
	volatile qreal mag = sqrt(unitAxis.x*unitAxis.x+unitAxis.y*unitAxis.y+unitAxis.z*unitAxis.z);
	unitAxis.x = unitAxis.x/mag;
	unitAxis.y = unitAxis.y/mag;
	unitAxis.z = unitAxis.z/mag;
	volatile qreal alphaReal = cos(g->angle/2.0);
	volatile qreal alphaImag = - sin(g->angle/2.0)*unitAxis.z;
	volatile qreal betaReal = sin(g->angle/2.0)*unitAxis.y;
	volatile qreal betaImag = - sin(g->angle/2.0)*unitAxis.x;

	for (thisTask=0; thisTask<numTasks; thisTask++) {
		enclose_offset = thisTask+offset+cacheid*numTasks;
		controlBit = extractBit (controlQubit, enclose_offset+chunkId*chunkSize);
		if (controlBit){
			// store current state vector values in temp variables
			stateRealUp = stateVecRealUp[enclose_offset];
			stateImagUp = stateVecImagUp[enclose_offset];

			stateRealLo = stateVecRealLo[enclose_offset];
			stateImagLo = stateVecImagLo[enclose_offset];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			stateVecRealUp[enclose_offset] = alphaReal*stateRealUp - alphaImag*stateImagUp 
				- betaReal*stateRealLo - betaImag*stateImagLo;
			stateVecImagUp[enclose_offset] = alphaReal*stateImagUp + alphaImag*stateRealUp 
				- betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			stateVecRealLo[enclose_offset] = betaReal*stateRealUp - betaImag*stateImagUp 
				+ alphaReal*stateRealLo + alphaImag*stateImagLo;
			stateVecImagLo[enclose_offset] = betaReal*stateImagUp + betaImag*stateRealUp 
				+ alphaReal*stateImagLo - alphaImag*stateRealLo;
		}
	} 
} 

void tGateHalf(GateObject *g, const int64_t cacheid)
{       
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;

    const long long int chunkSize=qureg.numAmpsPerChunk;
	const int rankIsUpper=chunkIsUpper(qureg.chunkId,chunkSize,targetQubit);

    qreal *stateVecRealLo = qureg.stateVec.real;
    qreal *stateVecImagLo = qureg.stateVec.imag;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
	if (rankIsUpper) {
		offset = 0;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    qreal stateRealLo, stateImagLo;
    const qreal cosAngle = 1.0/sqrt(2);
    const qreal sinAngle = 1.0/sqrt(2);
	long long int enclose_offset;

	long long int index;
    for (index=0; index<numTasks; index++) {
		enclose_offset = index+offset+cacheid*numTasks;
		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		stateVecRealLo[enclose_offset] = cosAngle*stateRealLo - sinAngle*stateImagLo;
		stateVecImagLo[enclose_offset] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
    }
}
void sGateHalf(GateObject *g, const int64_t cacheid)
{       
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;

    const long long int chunkSize=qureg.numAmpsPerChunk;
	const int rankIsUpper=chunkIsUpper(qureg.chunkId,chunkSize,targetQubit);

    qreal *stateVecRealLo = qureg.stateVec.real;
    qreal *stateVecImagLo = qureg.stateVec.imag;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
	if (rankIsUpper) {
		offset = 0;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    qreal stateRealLo, stateImagLo;
    const qreal cosAngle = 0;
    const qreal sinAngle = 1;
	long long int enclose_offset;

	long long int index;
    for (index=0; index<numTasks; index++) {
		enclose_offset = index+offset+cacheid*numTasks;
		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		stateVecRealLo[enclose_offset] = cosAngle*stateRealLo - sinAngle*stateImagLo;
		stateVecImagLo[enclose_offset] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
    }
}

void pauliZHalf(GateObject *g, const int64_t cacheid)
{       
	Qureg qureg = g->q;
	int targetQubit = g->targetQubit;

    const long long int chunkSize=qureg.numAmpsPerChunk;
	const int rankIsUpper=chunkIsUpper(qureg.chunkId,chunkSize,targetQubit);

    qreal *stateVecRealLo = qureg.stateVec.real;
    qreal *stateVecImagLo = qureg.stateVec.imag;
	long long int offset = qureg.numAmpsPerChunk>>1LL;
	if (rankIsUpper) {
		offset = 0;
		stateVecRealLo = qureg.pairStateVec.real;
		stateVecImagLo = qureg.pairStateVec.imag;
	}

    qreal stateRealLo, stateImagLo;
    const qreal cosAngle = -1;
    const qreal sinAngle = 0;
	long long int enclose_offset;

	long long int index;
    for (index=0; index<numTasks; index++) {
		enclose_offset = index+offset+cacheid*numTasks;
		stateRealLo = stateVecRealLo[enclose_offset];
		stateImagLo = stateVecImagLo[enclose_offset];

		stateVecRealLo[enclose_offset] = cosAngle*stateRealLo - sinAngle*stateImagLo;
		stateVecImagLo[enclose_offset] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
    }
}

#ifdef __cplusplus
}
#endif

#endif //QUEST_CACHE_H
