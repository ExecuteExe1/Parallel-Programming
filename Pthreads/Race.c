#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include <unistd.h>
  

int mails=0;
void* routine(){
    for(int i=0;i<100000;i++){
        mails++;
        //read mails
        //increment the value
        //then you write mails
        //this is fine for 1 thread,now we use 2 threads
    }
 
}

int main(int argc,char* argv[]){
 pthread_t p1,p2;
 if(pthread_create(&p1,NULL,&routine,NULL)!=0){
     return 1;
 }
 if(pthread_create(&p2,NULL,&routine,NULL)!=0){
     return 2;
 }
 if(pthread_join(p1,NULL)!=0){
     return 3;
  }
  if(pthread_join(p2,NULL)!=0){
     return 4;
  }
 printf("Number of mails: %d\n",mails);
 return 0;
}

// we do not get right results if we increase the range from 100 to 100000 race condition!

// f.e       |t1|t2
// read mails|23|24 
// increment |23|24
//write mails|23|25

// f.e       |t1|t2
// read mails|23|23  <- what if t1 and t2 read 23 at the same time 
// increment |23|2
//write mails|24|2
//             ^
//             |
//       t1 continues till the write while t2 is paused on 23
//           |t1|t2
// read mails|24|23  <---t1 continues,we even get to 30 or something
// increment |24|
//write mails|25|

//           |t1|t2
// read mails|29|23  <---t2 suddenly awakes and starts
// increment |29|30
//write mails|30|24   <-- we went from 30 back to 24 ??? No correct results
