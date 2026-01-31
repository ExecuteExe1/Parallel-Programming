
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
 

int main(int argc,char* argv[]){
   int pid=fork();
    int x=2;
    if(pid == -1){
        return 1;
    }
     if(pid==0){
         x++;
     }
    sleep(2);
    printf("Value of x: %d\n ",x);
    if(pid!= 0){
        wait(NULL);
     }
    return 0;
 
 }



