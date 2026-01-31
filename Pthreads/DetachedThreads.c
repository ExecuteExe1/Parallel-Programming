
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#define THREAD_NUM 2

 void* routine(void* args){

 sleep(1);
 printf("Finished execution\n");

 }



int main(int argc, char* argv[]){

 pthread_t th[THREAD_NUM];
 pthread_attr_t detachedThreads;
 pthread_attr_init(&detachedThreads);
 pthread_attr_setdetachstate(&detachedThreads,PTHREAD_CREATE_DETACHED);
  for(int i=0;i<THREAD_NUM;i++){
      if(pthread_create(&th[THREAD_NUM],&detachedThreads,&routine,NULL)!=0){
         perror("Failed to create thread");
       }
       pthread_detach(th[i]);
  }
  pthread_attr_destroy(&detachedThreads);
   
   

//    for(int i=0;i<THREAD_NUM;i++){
 //    if(pthread_join(th[THREAD_NUM],NULL)!=0){
 //        perror("Failed to join thread");
 //    }
  
 //   }

   pthread_exit(0);  
 }