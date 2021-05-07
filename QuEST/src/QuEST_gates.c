#include "QuEST_gates.h"
#include "QuEST_internal.h"
#include "QuEST_cache.h"
#include <stdlib.h>
#include <string.h>


// qureg ext handling


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
    return ext;
}


void calcOutcome_firstrun(Qureg *q) {
    if (q->ext->firstrun)
        q->ext->firstrun = 0;
    else return;
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
    register_gate_object(q, g);
}
