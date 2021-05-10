#include "QuEST_gates.h"
#include "QuEST_internal.h"
#include "QuEST_cache.h"
#include "QuEST_distribute_cache.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>


// qureg ext handling
// write the statevc into a file
static void write_file_f(Qureg *qureg) {
	FILE *f=fopen("midfile_f.dat", "w");
	long long int total_amp = qureg->numAmpsPerChunk;
	for (long long int i=0; i<total_amp; ++i)
		fprintf(f, "[%ld]\t%12.6f : %12.6f\n", i, qureg->pairStateVec.real[i], qureg->pairStateVec.imag[i]);
}
static void write_file_l(Qureg *qureg) {
	FILE *f=fopen("midfile_l.dat", "w");
	long long int total_amp = qureg->numAmpsPerChunk;
	for (long long int i=0; i<total_amp; ++i)
		fprintf(f, "[%ld]\t%12.6f : %12.6f\n", i, qureg->pairStateVec.real[i], qureg->pairStateVec.imag[i]);
}
static inline void clear_gate(QuregExt *ext) {
	GateObject *cur = ext->head;
	while(cur) {
		GateObject *tmp = cur;
		cur = cur->next;
		free(tmp);
	}
	ext->head = ext->last = NULL;
}

QuregExt* init_quregext(Qureg *q) {
    // normal init
    QuregExt *ext = malloc(sizeof(QuregExt));
    ext->firstrun = 1;
    ext->head = ext->last = NULL;
    // calculate num of local qubits
    long long int numtasks = q->numAmpsPerChunk;
    int maxbits = 0;
    while ((1LL << (maxbits + 1)) <= numtasks) maxbits++;
    ext->chunkWidth = maxbits;
	ext->exchanged = 0;
    return ext;
}

// judge if require_num gates' targetQubit are the same
int same_targetQubit(GateObject *cur, int require_num) {
	GateObject *tmp = cur->next;
	int targetQubit = cur->targetQubit;
	while (tmp && require_num && tmp->targetQubit==targetQubit && tmp->func3) {
		require_num--;
		tmp = tmp->next;
	}
	return !require_num;
}

// exchange half of the memory in one process
void exchange_half_mem(Qureg *q, int targetQubit) {
	double t1 = MPI_Wtime();
    int rankIsUpper = chunkIsUpper(q->chunkId, q->numAmpsPerChunk, targetQubit);
	int pairRank = getChunkPairId(rankIsUpper, q->chunkId, q->numAmpsPerChunk, targetQubit);
	int TAG = 200;
	MPI_Status status;
	long long int maxMessageCount = (1LL<<28);
	long long int half_mem = q->numAmpsPerChunk >> 1;
	if (half_mem < maxMessageCount)
		maxMessageCount = half_mem;
	int numMessages = half_mem / maxMessageCount;
	int i;
	long long int send_offset = 0;
	long long int recv_offset = half_mem;
	if (rankIsUpper) {
		send_offset = half_mem; // send the lower part to the lower rank
		recv_offset = 0; // recieve the upper part from the lower rank
	}
	// exchange data
	for (i = 0; i < numMessages; ++i) {
        MPI_Sendrecv(&q->stateVec.real[send_offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &q->pairStateVec.real[recv_offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&q->stateVec.imag[send_offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &q->pairStateVec.imag[recv_offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
		send_offset += maxMessageCount;
		recv_offset += maxMessageCount;
	}
	double t2 = MPI_Wtime();
	printf("pid = %d, pairRank = %d, time = %12.6fs (transfer half)\n", q->chunkId, pairRank, t2 - t1);
}

// exchange half of the memory in one process back
void exchange_half_mem_back(Qureg *q, int targetQubit) {
	double t1 = MPI_Wtime();
    int rankIsUpper = chunkIsUpper(q->chunkId, q->numAmpsPerChunk, targetQubit);
	int pairRank = getChunkPairId(rankIsUpper, q->chunkId, q->numAmpsPerChunk, targetQubit);
	int TAG = 300;
	MPI_Status status;
	long long int maxMessageCount = (1LL<<28);
	long long int half_mem = q->numAmpsPerChunk >> 1;
	if (half_mem < maxMessageCount)
		maxMessageCount = half_mem;
	int numMessages = half_mem / maxMessageCount;
	int i;
	long long int send_offset = half_mem;
	long long int recv_offset = 0;
	if (rankIsUpper) {
		send_offset = 0; // send the upper part back to the lower rank
		recv_offset = half_mem; // recieve the lower part from the lower rank
	}
	// exchange data
	for (i = 0; i < numMessages; ++i) {
        MPI_Sendrecv(&q->pairStateVec.real[send_offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &q->stateVec.real[recv_offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&q->pairStateVec.imag[send_offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &q->stateVec.imag[recv_offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
		send_offset += maxMessageCount;
		recv_offset += maxMessageCount;
	}
	double t2 = MPI_Wtime();
	printf("pid = %d, pairRank = %d, time = %12.6fs (transfer half back)\n", q->chunkId, pairRank, t2 - t1);
}

void calcOutcome_firstrun(Qureg *q) {
    if (!(q->ext->head))  // the list is empty
		return;
    // call handlers
    for (GateObject *cur = q->ext->head; cur; ) {
        if (cur->targetQubit < q->ext->chunkWidth && cur->func2) {
            GateObject *end = cur->next;
            // #define NO_GATE_FUSION
#ifndef NO_GATE_FUSION
            if (cur->targetQubit < CACHE_WIDTH) {
                while (end && end->func2 && end->targetQubit < CACHE_WIDTH)
                    end = end->next;
            }
            else {
                while (end && end->func2 && end->targetQubit == cur->targetQubit)
                    end = end->next;
            }
#endif //NO_GATE_FUSION
            #pragma omp parallel for schedule(static) default(shared)
            for (long long int cid = 0; cid < q->numAmpsPerChunk/CACHE_SIZE; cid++) {
                for (GateObject *cur2 = cur; cur2 != end; cur2 = cur2->next)
                    cur2->func2(cur2, cid);
            }
            cur = end;
#ifndef NO_MULTI_PROC_GATE_FUSION
        } else if (cur->targetQubit >= q->ext->chunkWidth && same_targetQubit(cur, 2)){
			int same_targetqubit = cur->targetQubit; // exchange half memory
			exchange_half_mem(q, same_targetqubit);
			// computes
			GateObject *end = cur->next;
			while (end && end->targetQubit == same_targetqubit && end->func3)
				end = end->next;
			long long int loop = q->numAmpsPerChunk / CACHE_SIZE;
			#pragma omp parallel for schedule(static) default(shared)
			for (long long int cid = 0; cid < loop; cid++) {
				for (GateObject *cur2 = cur; cur2 != end; cur2 = cur2->next)
					cur2->func3(cur2, cid);
			}
			exchange_half_mem_back(q, same_targetqubit);
			cur = end;
#endif // NO_MULTI_PROC_GATE_FUSION
		} else {
            cur->func(cur);
            cur = cur->next;
        }
    }
    // calc prob for once
    if (!q->isDensityMatrix) {
        long long int numtasks = q->numAmpsPerChunk;
        qreal *svreal = q->stateVec.real;
        qreal *svimag = q->stateVec.imag;
        qreal *output = q->ext->proboutcome;
        int maxbits = q->ext->chunkWidth;
        memset(output, 0, (maxbits + 1) * sizeof(qreal));
        #pragma omp parallel for schedule(static) default(shared) \
                reduction(+:output[:maxbits+1])
        for (long long int i = 0; i < numtasks; i++) {
            qreal val = svreal[i]*svreal[i] + svimag[i]*svimag[i];
            for (int b = 0; b < maxbits; b++) {
                if (!((i>>b) & 1)) {
                    output[b] += val;
                }
            }
            output[maxbits] += val;
        }
    }
	clear_gate(q->ext);
}


void destroy_quregext(QuregExt *ext) {
    GateObject *cur = ext->head;
    while (cur) {
        GateObject *tmp = cur;
        cur = cur->next;
        free(tmp);
    }
    free(ext);
}


// utils


static inline GateObject* create_gate_object() {
    return malloc(sizeof(GateObject));
}


static inline void register_gate_object(Qureg *q, GateObject *g) {
    // setup qureg
    g->q = *q;
    g->chunkOffset = q->chunkId * q->numAmpsPerChunk;
    g->svreal = q->stateVec.real;
    g->svimag = q->stateVec.imag;
    // setup linked list
    g->next = NULL;
    if (q->ext->last) {
        q->ext->last->next = g;
        q->ext->last = g;
    } else {
        q->ext->head = q->ext->last = g;
    }
}


// gates


static void handler_hadamard(GateObject *g) {
    statevec_hadamard(g->q, g->targetQubit);
}

void register_hadamard(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_hadamard;
    g->func2 = hadamardCache;
	g->func3 = hadamardHalf;
    register_gate_object(q, g);
}


static void handler_tGate(GateObject *g) {
    statevec_tGate(g->q, g->targetQubit);
}

void register_tGate(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_tGate;
    g->func2 = tGateCache;
	g->func3 = tGateHalf;
    register_gate_object(q, g);
}


static void handler_sGate(GateObject *g) {
    statevec_sGate(g->q, g->targetQubit);
}

void register_sGate(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_sGate;
    g->func2 = sGateCache;
    g->func3 = sGateHalf;
    register_gate_object(q, g);
}


static void handler_pauliX(GateObject *g) {
    statevec_pauliX(g->q, g->targetQubit);
}

void register_pauliX(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_pauliX;
    g->func2 = pauliXCache;
    g->func3 = pauliXHalf;
    register_gate_object(q, g);
}


static void handler_pauliY(GateObject *g) {
    statevec_pauliY(g->q, g->targetQubit);
}

void register_pauliY(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_pauliY;
    g->func2 = pauliYCache;
    g->func3 = pauliYHalf;
    register_gate_object(q, g);
}


static void handler_pauliZ(GateObject *g) {
    statevec_pauliZ(g->q, g->targetQubit);
}

void register_pauliZ(Qureg *q, int target) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->func = handler_pauliZ;
    g->func2 = pauliZCache;
    g->func3 = pauliZHalf;
    register_gate_object(q, g);
}


static void handler_rotateX(GateObject *g) {
    statevec_rotateX(g->q, g->targetQubit, g->angle);
}

void register_rotateX(Qureg *q, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_rotateX;
    g->func2 = rotateXCache;
    g->func3 = rotateXHalf;
    register_gate_object(q, g);
}


static void handler_rotateY(GateObject *g) {
    statevec_rotateY(g->q, g->targetQubit, g->angle);
}

void register_rotateY(Qureg *q, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_rotateY;
    g->func2 = rotateYCache;
    g->func3 = rotateYHalf;
    register_gate_object(q, g);
}


static void handler_rotateZ(GateObject *g) {
    statevec_rotateZ(g->q, g->targetQubit, g->angle);
}

void register_rotateZ(Qureg *q, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_rotateZ;
    g->func2 = rotateZCache;
    g->func3 = rotateZHalf;
    register_gate_object(q, g);
}


static void handler_controlledNot(GateObject *g) {
    statevec_controlledNot(g->q, g->controlQubit, g->targetQubit);
}

void register_controlledNot(Qureg *q, int control, int target) {
    GateObject *g = create_gate_object();
    g->controlQubit = control;
    g->targetQubit = target;
    g->func = handler_controlledNot;
    g->func2 = controlledNotCache;
    g->func3 = controlledNotHalf;
    register_gate_object(q, g);
}


static void handler_controlledPauliY(GateObject *g) {
    statevec_controlledPauliY(g->q, g->controlQubit, g->targetQubit);
}

void register_controlledPauliY(Qureg *q, int control, qreal target) {
    GateObject *g = create_gate_object();
    g->controlQubit = control;
    g->targetQubit = target;
    g->func = handler_controlledPauliY;
    g->func2 = controlledPauliYCache;
    g->func3 = controlledPauliYHalf;
    register_gate_object(q, g);
}


static void handler_controlledRotateX(GateObject *g) {
    statevec_controlledRotateX(g->q, g->controlQubit, g->targetQubit, g->angle);
}

void register_controlledRotateX(Qureg *q, int control, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->controlQubit = control;
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_controlledRotateX;
    g->func2 = controlledRotateXCache;
    g->func3 = controlledRotateXHalf;
    register_gate_object(q, g);
}


static void handler_controlledRotateY(GateObject *g) {
    statevec_controlledRotateY(g->q, g->controlQubit, g->targetQubit, g->angle);
}

void register_controlledRotateY(Qureg *q, int control, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->controlQubit = control;
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_controlledRotateY;
    g->func2 = controlledRotateYCache;
    g->func3 = controlledRotateYHalf;
    register_gate_object(q, g);
}


static void handler_controlledRotateZ(GateObject *g) {
    statevec_controlledRotateZ(g->q, g->controlQubit, g->targetQubit, g->angle);
}

void register_controlledRotateZ(Qureg *q, int control, int target, qreal angle) {
    GateObject *g = create_gate_object();
    g->controlQubit = control;
    g->targetQubit = target;
    g->angle = angle;
    g->func = handler_controlledRotateZ;
    g->func2 = controlledRotateZCache;
    g->func3 = controlledRotateZHalf;
    register_gate_object(q, g);
}
