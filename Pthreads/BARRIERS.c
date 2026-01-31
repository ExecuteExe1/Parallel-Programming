#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>

pthread_barrier_t barrier;
//barrier,makes multiple threads wait before letting MULTIPLE threads pass
void* routine(void* args){
  while(1){
    printf("Waiting at the barrier....\n");
   sleep(1);
   pthread_barrier_wait(&barrier);
   printf("We passed the barrier\n");
   sleep(1);
  }
}

int main(int argc,char* argv[]){

 pthread_t th[10];
 pthread_barrier_init(&barrier,NULL,3);
 int i;
 for(i=0;i<10;i++){
   if(pthread_create(&th[i],NULL,&routine,NULL)!=0){
     perror("Failed to create thread");
   }
 }
  for(i=0;i<10;i++){
   if(pthread_join(th[i],NULL)!=0){
     perror("Failed to join thread");
   }
 }
 pthread_barrier_destroy(&barrier);




return 0;
}