#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>


#include <signal.h>
#include <sys/wait.h>

int prepare(void);

int process_arglist(int count, char** arglist);

int finalize(void);

void handle_sigchld(int sig);
void reset_SIGCHILD(void);
void change_SIGCHILD(void);

void ignore_SIGINT(void);
void reset_SIGINT(void);

/* handle_sigchld taken from here: http://www.microhowto.info/howto/reap_zombie_processes_using_a_sigchld_handler.html this function will be called when our shell get SIGCHLD from a child process and call to waitpid if this child proccess if and only if terminated */

void handle_sigchld(int sig) {
       int saved_errno = errno;
       while (waitpid((pid_t)(-1), 0, WNOHANG) > 0) {}
       errno = saved_errno;
}
void change_SIGCHILD(void){
    /* code taken from here: http://www.microhowto.info/howto/reap_zombie_processes_using_a_sigchld_handler.html
     */
    
    struct sigaction sa;
    sa.sa_handler = &handle_sigchld;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_NOCLDSTOP;
    if (sigaction(SIGCHLD, &sa, 0) == -1) {
        printf("Error: %s\n", strerror(errno));
        exit(1);
    }
}


void reset_SIGCHILD(void){
    struct sigaction sa = {.sa_handler = SIG_DFL};
    if (sigaction(SIGCHLD, &sa, 0) == -1) {
        printf("Error: %s\n", strerror(errno));
        exit(1);
    }
}

    
void ignore_SIGINT(void){
    struct sigaction sa = {.sa_handler = SIG_IGN};
    if (sigaction(SIGINT, &sa, 0) == -1) {
        printf("Error: %s\n", strerror(errno));
        exit(1);
    }
}
    
void reset_SIGINT(void){
    struct sigaction sa = {.sa_handler = SIG_DFL};
    if (sigaction(SIGINT, &sa, 0) == -1) {
        printf("Error: %s\n", strerror(errno));
        exit(1);
    }
}



int prepare(void){
    change_SIGCHILD();
    ignore_SIGINT();
    
    return 0;
}


/*
 main idea:
 
 we have 4 distinguished cases. we will label each case by 0-3
 0: for rag command
 1: for command ended by &
 2: the >> command
 3: the pipe command
 

 
 the first step of process_arglist() would be to find the label of the command
 at this step we also modify some of the commands:
 
 for case 0 we dont channge the command
 for case 1, 2 and 3 we will raplace '&', '>>', '|' with NULL
 
 we use functions to change the signal handler as needed per case and per procces
 
 */

int process_arglist(int count, char** arglist){
    int i;
    int rc;
    int rc2;
    int label = 0;                     /* default case is regular command that we need to wait for */
    char** arglist2 = arglist;         /*just for getting clean compiling*/
    int pipefd[2];
    int fd;
    int status;
    
    if(strcmp("&", arglist[count-1])==0){
        label = 1;
        arglist[count-1] = arglist[count];
    }
    
    else{
        for(i=0;  i<count; i++){
            if(strcmp(">>", arglist[i])==0){
                label = 2;
                arglist[i] = arglist[count];
                break;
            }
            else if(strcmp("|", arglist[i])==0){
                /* the aim is to replace '|' with null, now arglist == part1+null+part2+null and arglist2 will point to part2+null*/
                label = 3;
                arglist[i] = arglist[count];
                arglist2 = arglist + i + 1;
                break;
            }
        }
    }
               
    if( label==3 ){         /*open a pipe if needed*/
        if(pipe(pipefd)==-1){
            printf("Error: %s\n", strerror(errno));
            exit(1);
        }
    }

    rc = fork();  /*  create new duplicated  proccess  */
    
    if(rc<0){   /* enter if fork failed */
        
        printf("Error: %s\n", strerror(errno));
        exit(1);
    }
    
    else if(rc==0){             /* will run this code only for the new proccess*/
        
        reset_SIGCHILD();
        reset_SIGINT();
        
        if( label ==0){
            
            execvp(arglist[0], arglist);
            printf("Error: %s\n", strerror(errno));
            exit(1);              /*this will execute only if execvp fails!*/
        }
        
        else if(label ==1){             /* & case, we need to ignore sigint */
            
            ignore_SIGINT();
            
            execvp(arglist[0], arglist);
            printf("Error: %s\n", strerror(errno));
            exit(1);
        }
        else if(label == 2){             /* >> case */
            
            fd = open(arglist[count-1], O_CREAT|O_WRONLY|O_APPEND, S_IRUSR|S_IWUSR);
            if(fd==-1) exit(1);
            
            if(dup2(fd, 1)==-1){        /*modify to write to file as defult*/
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            
            if(close(fd)==-1){  /*we only decrement the refrence count of the open file to 1*/
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            
            execvp(arglist[0], arglist);
            printf("Error: %s\n", strerror(errno));
            exit(1);
        }
        else{   /* child with the writing side of the pipe*/
            
            if( close(pipefd[0])==-1){            /*we dont need it*/
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            if(dup2(pipefd[1], 1)==-1){          /*duflt writing to the pipe*/
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            if( close(pipefd[1])==-1){          /*we only decrement the refrence count of the open file to 1*/
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            
            execvp(arglist[0], arglist);
            
            printf("Error: %s\n", strerror(errno)); /*in case execvp fails*/
            exit(1);
        }
                
    }
    
    else{   /* only parent will run this code after spliting from first fork call */
        
        if(label == 0 || label == 2){
            if(waitpid(rc, &status, 0)==-1){
                if(errno!=EINTR && errno!=ECHILD){
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
            } return 0;
        }
        /* in the case of label==1 we do noting so its either if or else if*/
        else if(label == 3){
            
            rc2 = fork();
            
            if(rc2<0){   /* enter if fork failed */
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            
            if(rc2==0){  /*2nd child of the pipe case, this child will read from the pipe*/
                
                reset_SIGCHILD();
                reset_SIGINT();
                
                if( close(pipefd[1]) ==-1){
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
                if( dup2(pipefd[0], 0) ==-1){
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
                if( close(pipefd[0]) ==-1){   /*as before we only decrement*/
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
                
            
                execvp(arglist2[0], arglist2);
                printf("Error: %s\n", strerror(errno));
                exit(1);
            }
            
            else{         /* we are inside perent in the case of pipe */
                
                if( close(pipefd[1])==-1){            /*we dont need it*/
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
                if( close(pipefd[0])==-1){            /*we dont need it*/
                    printf("Error: %s\n", strerror(errno));
                    exit(1);
                }
                
                if(waitpid(rc, &status, 0)==-1){
                    if(errno!=EINTR && errno!=ECHILD){
                        printf("Error: %s\n", strerror(errno));
                        exit(1);
                    }
                }
                
                if(waitpid(rc2, &status, 0)==-1){
                    if(errno!=EINTR && errno!=ECHILD){
                        printf("Error: %s\n", strerror(errno));
                        exit(1);
                    }
                } return 0;
            }
        }
        
    }
    return 0;
}

int finalize(void){
    return 0;
}


