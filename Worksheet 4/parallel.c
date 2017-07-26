#include "parallel.h"
#include  <math.h>

void init_parallel (
                    int iproc,
                    int jproc,
                    int imax,
                    int jmax,
                    int *myrank,
                    int *il,
                    int *ir,
                    int *jb,
                    int *jt,
                    int *rank_l,
                    int *rank_r,
                    int *rank_b,
                    int *rank_t,
                    int *omg_i,
                    int *omg_j,
                    int num_proc
                    ){

    // Get "domain size" for each subdomain
    int imax_sub = imax / iproc;
    int jmax_sub = jmax / jproc;

    // Get rank of this process
    MPI_Comm_rank(MPI_COMM_WORLD, myrank);
    
    // Get indices of each subdomain (refer to Figure 2 on worksheet 4) - This will provide subdomain index that corresponds to each taskID (rank)
    // Instead of getting the omg_i,j of all the subdomains and sending them to each rank, the corresponding subdomain indices are computed by each rank directly
    // Each Omega coordinates starts from 0 to be consistent with rank enumeration

    //Note that modulus in C follows the equation: x % y = x - floor(x/y)*y
    *omg_i = *myrank % iproc;
    
    //Note that integer division in C follows: x / y = floor(x/y) if x and y are integers
    *omg_j = *myrank / iproc;

    // Get left neighbour ranks 
    if (*omg_i > 0 && *omg_i < iproc){
        *rank_l = *myrank - 1;
    }
    else{    
        *rank_l = MPI_PROC_NULL;           
    }

    // Get right neighbour ranks
    if (*omg_i >= 0 && *omg_i < iproc - 1){ 
        *rank_r = *myrank + 1;
    }
    else{
        *rank_r = MPI_PROC_NULL;
    }

    // Get lower neighbour ranks
    if (*omg_j > 0 && *omg_j < jproc){
        *rank_b = *myrank - iproc;
    }
    else{
        *rank_b = MPI_PROC_NULL;
    }

    // Get upper neighbour ranks
    if (*omg_j >= 0 && *omg_j < jproc - 1){
        *rank_t = *myrank + iproc;
    }
    else{
        *rank_t = MPI_PROC_NULL;
    }
   
    // Get subdomain dimension for each omega (rank)
    *il = *omg_i * imax_sub + 1;
    *ir = *il + imax_sub - 1;
    *jb = *omg_j * jmax_sub + 1;
    *jt = *jb + jmax_sub - 1;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Please note MPI_PROC_NULL will be printed as negative numbers. Condition statements used in Communication.c could be conisdered but it gets too messy here...
    printf("TaskID=%d treates subdomain Omega(%d, %d) with size =[%d, %d]x[%d, %d] || Neighbour TaskID: Left=%d Right=%d Down=%d Up=%d\n", *myrank, *omg_i, *omg_j, *il, *ir, 
            *jb, *jt, *rank_l, *rank_r, *rank_b, *rank_t);
    printf("\n");

    MPI_Barrier(MPI_COMM_WORLD);

}

/* Function void pressure_comm */
/* this function will be called by SOR before computing residual */
void pressure_comm(double **P,
		           int  il,
		           int ir,
 		           int jb,
  		           int jt,
		           int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
       		       double *bufRecv
                   //MPI_Status *status -not sure what to use this for
		           //int chunk  - not sure what this is meant to be doing?
                   ){

    //Here we could add further condition statements for MPI_NULL_PROC. 
    communicate(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, LEFT, VAR_P, bufSend, bufRecv);
    communicate(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, RIGHT, VAR_P, bufSend, bufRecv);
    communicate(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, UP, VAR_P, bufSend, bufRecv);
    communicate(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, DOWN, VAR_P, bufSend, bufRecv);

}


void uv_comm(
            double **U,
            double **V,
            int il, 
            int ir,
            int jb, 
            int jt,
            int rank_l, 
            int rank_r,
            int rank_b, 
            int rank_t,
            double *bufSend,
            double *bufRecv 
           ){ 
     
     communicate(U, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, RIGHT, VAR_U, bufSend, bufRecv);
     communicate(V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, RIGHT, VAR_V, bufSend, bufRecv);
 
     communicate(U, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, LEFT, VAR_U, bufSend, bufRecv);
     communicate(V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, LEFT, VAR_V, bufSend, bufRecv);
     
     communicate(U, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, UP, VAR_U, bufSend, bufRecv);
     communicate(V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, UP, VAR_V, bufSend, bufRecv);
 
     communicate(U, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, DOWN, VAR_U, bufSend, bufRecv);
     communicate(V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, DOWN, VAR_V, bufSend, bufRecv);
     
}                  


/* produces a stderr text output  */
void Program_Message(char *txt){

   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);

}


/* produces a stderr textoutput and synchronize all processes */
void Programm_Sync(char *txt){

   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);

}


/* all processes will produce a text output, be synchronized and finished */
void Programm_Stop(char *txt){

   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);

}


//CHECK_NP checks if the np value given is appropriate. If np is an odd number the program will stop (this is the limitation of our current solver)
int check_np (
               int *iproc,
               int *jproc,
               int num_proc,
               int myrank
             ){

    //Check if number of processes is even (odd number won't work for this worksheet)
     if (num_proc % 2 != 0){

        if (num_proc == 1){

            *iproc = 1;
            *jproc = 1;

            if (myrank == 0){

            printf("\n");
            printf("Number of processes = %d\n", num_proc);
            printf("\n");
            printf("Program will solve the problem sequentially (iproc = %d, jproc = %d)\n", *iproc, *jproc);
            printf("\n");
    
            }
    
        }

        else{

            if (myrank == 0){

            printf("\n");
            printf("Inappropriate number of processes (please use even number of processes)\n");
            printf("\n");
            printf("Exiting program...\n");
            printf("\n");

            }

            MPI_Finalize();
            exit(EXIT_FAILURE);

        }

    }

    if(*iproc * *jproc != num_proc){

        if (myrank == 0){

        printf("\n");
        printf("Incompatible number of processes: iproc=%d, jproc=%d, np=%d (np must be iproc * jproc)\n", *iproc, *jproc, num_proc);
        printf("\n");
    
        }

        correct_ijproc(iproc, jproc, num_proc);

        if (myrank == 0){

        printf("iproc and jproc redefined: iproc = %d, jproc = %d for np = %d\n", *iproc, *jproc, num_proc);
        printf("\n");
        
        }

    }

    else{
    
    return 0;
    
    }

    return 0;

}

//CORRECT_IJPROC redefines the iproc and jproc values based on the given np value (if this does not match with the original iproc and jproc values)
void correct_ijproc (
                     int *iproc,
                     int *jproc,
                     int num_proc
                    ){

    double sqroot;

    sqroot = sqrt(num_proc);

    if ((sqroot - floor(sqroot)) == 0){

        *iproc = (int) sqroot;
        *jproc = (int) sqroot;

    }

    else{

        *iproc = (int) num_proc / 2;
        *jproc = (int) num_proc / *iproc;

    }

}


//send_buffer sends the data in MATRIX with boundary il, ir, jb, jt to buffer space. Note that **matrix corresponds to each subdomain
void send2buffer(
                 double **matrix, 
                 double *buffer, 
                 int il, 
                 int ir, 
                 int jb, 
                 int jt
                ){

    int i, j;
    int k = 0;

    for(i = il; i <= ir; i++){

        for(j = jb; j <= jt; j++){

            buffer[k++] = matrix[i][j];

        }

    }

}

//recv_buffer retrieves data saved in the buffer space to the corresponding rank (name recv2buffer is just for consistency. Ignore the grammar!)
//Note that **matrix corresponds to each subdomain
void recv2buffer(
                 double **matrix, 
                 double *buffer, 
                 int il, 
                 int ir, 
                 int jb, 
                 int jt
                ){

    int i, j;
    int k = 0;

    for(i = il; i <= ir; i++){

        for(j = jb; j <= jt; j++){

            matrix[i][j] = buffer[k++];

        }

    }

}    


//COMMUNICATE exchanges data between two ranks using buffer space
//Some idea referred from many great online MPI tutorials including : https://computing.llnl.gov/tutorials/mpi/
void communicate(
                 double **matrix,
                 int il, 
                 int ir, 
                 int jb, 
                 int jt,
                 int rank_l, 
                 int rank_r, 
                 int rank_b, 
                 int rank_t,
                 int direction,
                 int variable,
                 double *bufSend,
                 double *bufRecv
                ){

    int num_elem;                                                       //number of elements to be sent/received
    int destination, source;                                            //rank the information is being sent from/to
    int send_i_min, send_i_max, send_j_min, send_j_max;                 //"range" of subdomain indices for data to be sent
    int receive_i_min, receive_i_max, receive_j_min, receive_j_max;     //"range" of subdomain indices to which the data is to be retrieved
    int count_received;                                                 //number of elements being received (only acts as safety switch here) 

    MPI_Status status;
    
    switch(direction){

        case RIGHT:

            num_elem = jt - jb + 3;

            destination = rank_r;
            source = rank_l;

            send_i_min = ir;
            send_i_max = ir;
            send_j_min = jb;
            send_j_max = jt;

            receive_i_min = il - 1;
            receive_i_max = il - 1;
            receive_j_min = jb;
            receive_j_max = jt;

            switch(variable){

                case VAR_P:

                    send_i_min = ir;
                    send_i_max = ir;
                    send_j_min = jb;
                    send_j_max = jt;
    
                    receive_i_min = il - 1;
                    receive_i_max = il - 1;
                    receive_j_min = jb;
                    receive_j_max = jt;
        
                    break;

                case VAR_V:

                    send_i_min = ir;
                    send_i_max = ir;
                    send_j_min = jb - 1;
                    send_j_max = jt;

                    receive_i_min = il - 1;
                    receive_i_max = il - 1;
                    receive_j_min = jb - 1;
                    receive_j_max = jt;

                    break;

                case VAR_U:

                    send_i_min = ir;
                    send_i_max = ir;
                    send_j_min = jb;
                    send_j_max = jt;

                    receive_i_min = il - 2;
                    receive_i_max = il - 2;
                    receive_j_min = jb;
                    receive_j_max = jt;

                    break;

            }

            break;

        case LEFT:

            num_elem = jt - jb + 2;

            destination = rank_l;
            source = rank_r;

            send_i_min = il;
            send_i_max = il;
            send_j_min = jb;
            send_j_max = jt;
 
            receive_i_min = ir + 1;
            receive_i_max = ir + 1;
            receive_j_min = jb;
            receive_j_max = jt;                                                     

            switch(variable){
                     
                case VAR_P:
                     
                     send_i_min = il;
                     send_i_max = il;
                     send_j_min = jb;
                     send_j_max = jt;
                     
                     receive_i_min = ir + 1;
                     receive_i_max = ir + 1;
                     receive_j_min = jb;
                     receive_j_max = jt;
                     
                     break;

                case VAR_U:

                    send_i_min = il - 1;
                    send_i_max = il - 1;
                    send_j_min = jb;
                    send_j_max = jt;
        
                    receive_i_min = ir + 1;
                    receive_i_max = ir + 1;
                    receive_j_min = jb;
                    receive_j_max = jt;

                    break;

                case VAR_V:

                    send_i_min = il;
                    send_i_max = il;
                    send_j_min = jb - 1;
                    send_j_max = jt;
    
                    receive_i_min = ir + 1;
                    receive_i_max = ir + 1;
                    receive_j_min = jb - 1;
                    receive_j_max = jt;

                    break;

            }
 
            break;

        case UP:

            num_elem = ir - il + 2;

            destination = rank_t;
            source = rank_b;

            send_i_min = il;
            send_i_max = ir;
            send_j_min = jt;
            send_j_max = jt;

            receive_i_min = il;
            receive_i_max = ir;
            receive_j_min = jb - 1;
            receive_j_max = jb - 1;

            switch(variable) {

                case VAR_P:

                    send_i_min = il;
                    send_i_max = ir;
                    send_j_min = jt;
                    send_j_max = jt;

                    receive_i_min = il;
                    receive_i_max = ir;
                    receive_j_min = jb - 1;
                    receive_j_max = jb - 1;
    
                    break;

                case VAR_U:

                    send_i_min = il - 1;
                    send_i_max = ir;
                    send_j_min = jt;
                    send_j_max = jt;

                    receive_i_min = il - 1;
                    receive_i_max = ir;
                    receive_j_min = jb - 1;
                    receive_j_max = jb - 1;

                    break;

                case VAR_V:

                    send_i_min = il;
                    send_i_max = ir;
                    send_j_min = jt;
                    send_j_max = jt;

                    receive_i_min = il;
                    receive_i_max = ir;
                    receive_j_min = jb - 2;
                    receive_j_max = jb - 2;

                    break;

            }

            break;

        case DOWN:

            num_elem = ir - il + 2;

            destination = rank_b;
            source = rank_t;

            send_i_min = il;
            send_i_max = ir;
            send_j_min = jb;
            send_j_max = jb;

            receive_i_min = il;
            receive_i_max = ir;
            receive_j_min = jt + 1;
            receive_j_max = jt + 1;

        switch (variable){

            case VAR_P:

                 send_i_min = il;
                 send_i_max = ir;
                 send_j_min = jb;
                 send_j_max = jb;

                 receive_i_min = il;
                 receive_i_max = ir;
                 receive_j_min = jt + 1;
                 receive_j_max = jt + 1;

            break;
 
            case VAR_U:
 
                 send_i_min = il - 1;
                 send_i_max = ir;
                 send_j_min = jb;
                 send_j_max = jb;

                 receive_i_min = il - 1;
                 receive_i_max = ir;
                 receive_j_min = jt + 1;
                 receive_j_max = jt + 1;
 
             break;

             case VAR_V:

                  send_i_min = il;
                  send_i_max = ir;
                  send_j_min = jb;
                  send_j_max = jb;

                  receive_i_min = il;
                  receive_i_max = ir;
                  receive_j_min = jt + 1;
                  receive_j_max = jt + 1;

             break;
 
             }

            break;

    }

    send2buffer(matrix, bufSend, send_i_min, send_i_max, send_j_min, send_j_max);

    MPI_Sendrecv(bufSend, num_elem, MPI_DOUBLE, destination, 0, bufRecv, num_elem, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);

    MPI_Get_count(&status, MPI_DOUBLE, &count_received);

    if(count_received > 0){

        recv2buffer(matrix, bufRecv, receive_i_min, receive_i_max, receive_j_min, receive_j_max);    
    }
    
}


