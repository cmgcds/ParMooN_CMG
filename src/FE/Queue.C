#include <Queue.h>
void enqueue(int k, queue * qs){
    
  node * newnode = new node[1];
  newnode->next = NULL;
  newnode-> id = k;
  
  node *temp = qs->bottom;;
  if(temp){
    temp->next = newnode;
    qs->bottom = newnode;
  }
  else{
    qs->top = newnode;
    qs->bottom = newnode;
  }
}

int dequeue(queue *qs){
  int key=-1;
  node * temp=qs->top;
 
  if(temp){
    key = temp->id;
    if(temp->next)
      qs->top = temp->next;
    else{
      qs->top = NULL;
      qs->bottom = NULL;
    }
    
    temp = NULL;
    delete temp;        
  }  
  return key;
}

queue* createQ(){
  
  queue* k = new queue[1];
  
  k->top = NULL;
  k->bottom = NULL;
  
  return k;
  
}

node* createNode(int key){
  node * newnode = new node[1];
  newnode->id = key;
  newnode->next = NULL;
}