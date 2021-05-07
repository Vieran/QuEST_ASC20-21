#ifndef QUEST_GATES_H
#define QUEST_GATES_H
#include "QuEST.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct GateObject {
    // qureg
    Qureg q;
    int64_t chunkOffset;
    qreal *svreal, *svimag;
    // args
    int controlQubit, targetQubit;
    qreal angle;
    // others
    void (*func)(struct GateObject*);
    void (*func2)(struct GateObject*, const int64_t);
    struct GateObject* next;
} GateObject;


typedef struct QuregExt {
    int firstrun;
    int chunkWidth;
    qreal proboutcome[65];
    struct GateObject *head, *last;
} QuregExt;


QuregExt* init_quregext(Qureg *q);

void calcOutcome_firstrun(Qureg *q);

void destroy_quregext(QuregExt *qe);


void register_hadamard(Qureg *q, int target);

void register_tGate(Qureg *q, int target);

void register_sGate(Qureg *q, int target);

void register_pauliX(Qureg *q, int target);

void register_pauliY(Qureg *q, int target);

void register_pauliZ(Qureg *q, int target);

void register_rotateX(Qureg *q, int target, qreal angle);

void register_rotateY(Qureg *q, int target, qreal angle);

void register_rotateZ(Qureg *q, int target, qreal angle);

void register_controlledNot(Qureg *q, int control, int target);

void register_controlledPauliY(Qureg *q, int control, qreal target);

void register_controlledRotateX(Qureg *q, int control, int target, qreal angle);

void register_controlledRotateY(Qureg *q, int control, int target, qreal angle);

void register_controlledRotateZ(Qureg *q, int control, int target, qreal angle);


#ifdef __cplusplus
}
#endif

#endif //QUEST_GATES_H