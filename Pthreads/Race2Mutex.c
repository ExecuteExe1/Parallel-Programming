#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include <unistd.h>
  //works with more threads as well so we have a safe solution
int mails=0;
int lock=0;
pthread_mutex_t mutex;

void* routine(){
    for(int i=0;i<20000000;i++){
        pthread_mutex_lock(&mutex); //we do not care about waiting,
        mails++;                  
        lock=0;
        pthread_mutex_unlock(&mutex);
    }
 
}

int main(int argc,char* argv[]){
 pthread_t p1,p2,p3;
 pthread_mutex_init(&mutex,NULL);     //address of mutex and other attributes ,we need todestroy mutex at the end
 if(pthread_create(&p1,NULL,&routine,NULL)!=0){
     return 1;
 }
 if(pthread_create(&p2,NULL,&routine,NULL)!=0){
     return 2;
 }
 if(pthread_create(&p3,NULL,&routine,NULL)!=0){
     return 3;
 }
 if(pthread_join(p1,NULL)!=0){
     return 4;
  }
  if(pthread_join(p2,NULL)!=0){
     return 5;
  }
  if(pthread_join(p3,NULL)!=0){
     return 6;
  }
   pthread_mutex_destroy(&mutex); //destry memory
 printf("Number of mails: %d\n",mails);
 return 0;
}