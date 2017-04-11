

static long tube_1_init(aSubRecord *prec) {
  //initalization of asubrecord (possibly not necessary)


  printf("Initializing asub:\n");
  printf("number of elements in a:%i\n", prec->noa );
  printf("number of elements in b:%i\n", prec->nob );
  
  long i;
  double *a;
  a = (double *)prec->a;
  for (i=0; i<prec->noa; i++) {
    a[i]= (double)i/prec->noa;
    }

  long *b;
  b = (long *)prec->b;
  for (i=0; i<prec->nob; i++) {
    b[i]= i;
  }

  double *c;
  c = (double *)prec->c;
  for (i=0; i<prec->noc; i++) {
    c[i]= 0.25;
    printf("c initialized to %g\n", c[i]);
  }

  return 0;

}

static long tube_1_routine(aSubRecord *prec) {
  //get new data from digitizer, since last call
  //each tube (digitizer input) needs:
  //        array for energy
  //        array for time
  //        scaler for number of counts

  long i;
  double *Energy;
  double sum=0;
  Energy = (double *)prec->a;
  for (i=0; i<prec->noa; i++) {
    printf(" Energy[%lu]", i);
    printf( ":%g\n", Energy[i]);
    Energy[i]= Energy[i]+(double)i/prec->noa;
  }

  double *Time;
  Time = (double *)prec->b;
  for (i=0; i<prec->nob; i++) {
    printf(" Time[%lu]", i);
    printf( ":%lu\n", Time[i]);
    Time[i]= Time[i]*i;
  }
 
  double *c;
  scalar = (double *)prec->c;
  for (i=0; i<prec->noc; i++) {
    printf( "Value of scalar is: %g\n", c[i]);
    scalar[i]= scalar[i]*2;
  }

  //clear digitizer?
  return(0);
};


epicsRegisterFunction(tube_1_init);
epicsRegisterFunction(tube_1_routine);
