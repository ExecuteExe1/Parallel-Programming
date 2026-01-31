#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

 //getting the return value INTO our function

void* roll_dice(){
  int value=(rand()% 6)+ 1 ;
  int* result=malloc(sizeof(int));
 // printf("%d\n",value);
   printf("Thread result: %p\n", result);
   return (void*) result;
   }

int main(int argc,char* argv[]){
     int* res;
     pthread_t th;
     srand(time(NULL));
     if(pthread_create(&th,NULL,&roll_dice,NULL)!=0){
         return 1;
      }
     if(pthread_join(th,(void**) &res)!=0){
       return 2;
     }
     printf("Main res: %p\n",res);
     printf("Result:%d\n",*res);
      free(res);
     return 0;  
 }

