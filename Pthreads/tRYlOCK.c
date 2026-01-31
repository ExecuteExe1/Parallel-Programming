#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

pthread_mutex_t mutex;

void* routine(void* arg){
  if(pthread_mutex_trylock(&mutex)==0){  //not guranteed,lock is gonna wait for mutex to be unlocked,unlocked will NOT wait
   printf("Got lock\n");
   sleep(1);
   pthread_mutex_unlock(&mutex);
  }
  else{
    printf("Did not get lock\n");
  }
 
}


int main(int argc,char* argv[]){
 pthread_t th[4];
 pthread_mutex_init(&mutex,NULL);
  for(int i=0;i<4;i++){
     if(pthread_create(&th[i],NULL,&routine,NULL)!=0){
     perror("Error creating a thread");
 }

}
   for(int i=0;i<4;i++){
     if(pthread_join(th[i],NULL)!=0){
        perror("Error at joining thread");
    }
   }
 pthread_mutex_destroy(&mutex);
 return 0;
}