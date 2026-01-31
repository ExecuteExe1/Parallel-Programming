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
    for(int i=0;i<10000000;i++){
        pthread_mutex_lock(&mutex); //we do not care about waiting,
        mails++;                  
        lock=0;
        pthread_mutex_unlock(&mutex);
    }
 
}

int main(int argc,char* argv[]){
 pthread_t th[4];
 pthread_mutex_init(&mutex,NULL);   
   for(int i=0;i<4;i++)
   {

     if(pthread_create(&th[i],NULL,&routine,NULL)!=0){
        perror("Failed to create Thread");
     return 1;
    }
    printf("Thread %d has started\n",i);
    if(pthread_join(th[i],NULL)!=0){
     return 2;
     }
       printf("Thread %d has finished its execution\n",i);
   }
//for(int i=0;i<4;i++){
//
 //   if(pthread_join(th[i],NULL)!=0){
//     return 2;
//     }
 //      printf("Thread %d has finished its execution\n",i);
//}
// this ,if our thread_joint is in another loop things become parallel! use // to make the joint part comment and remove //these to try this code
  


   pthread_mutex_destroy(&mutex); 
 printf("Number of mails: %d\n",mails);
 return 0;


}

//we have no threads running at the same time........we want things to run in parallel

