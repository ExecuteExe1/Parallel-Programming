#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
  
   int x=2;
void *routine(){
   
 x++;
 sleep(2);
  printf("Value of x: %d/n",x);

}
void *routine2(){
    
     sleep(2);
   printf("Value of x: %d/n",x);

}
 void *routine3(){
     x=2;
   printf("Value of x: %d/n",x);

}
 

int main(int argc, char* argv[]){

 pthread_t t1,t2,t3;
 if( pthread_create(&t1,NULL,&routine,NULL)){
      return 1;
       }
  if( pthread_create(&t2,NULL,&routine2,NULL)){
      return 2;
     }
    if( pthread_create(&t3,NULL,&routine3,NULL)){
      return 3;
       } 

  if(pthread_join(t1,NULL)){
      return 4;
    }
    if(pthread_join(t2,NULL)){
      return 5;
    }
   if(pthread_join(t3,NULL)){
      return 6;
    }

  return 0;   //threads share memories,the 1st thread made x=3 so the second slept and then it's x=3 too
 }