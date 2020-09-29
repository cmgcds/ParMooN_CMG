
#ifndef __Queue__
#define __Queue__
#include <stdlib.h>

struct node_n{
    int id;
    struct node_n * next;
};

typedef struct node_n node;

struct q_ravi{
  
  node * top;
  node * bottom;
  
} ;

typedef q_ravi queue;

void enqueue(int k, queue * qs);
int dequeue(queue *qs);
queue* createQ();
node* createNode(int key);
#endif
