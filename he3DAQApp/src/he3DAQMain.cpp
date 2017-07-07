/*
 *
 *   If anyone has to maintain this code, I am truly and deeply sorry :(
 *
 */

#include <stddef.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <cstring>

#include <dirent.h>
#include <errno.h>

#include <subRecord.h>
#include <aSubRecord.h>
#include "stringoutRecord.h"

#include "epicsExit.h"
#include "epicsThread.h"
#include "iocsh.h"
#include "epicsMessageQueue.h"

#include <registryFunction.h>
#include <epicsExport.h>

//root
#include "TFile.h"
#include "TTree.h" 


epicsMessageQueue *mq1;
epicsThreadId myThread=0;
epicsMutex *plock;
int message;
double *Data;
int runNum=0;


struct timeval timeMark;

//number of triggers before the data is saved to file
int NTRIGS=5080;
int m_numTrigs[8]={0};
bool isRunning=false;


int numSwTrigs=0;
/*
//vectors to hold the data before it gets saved to file
std::vector<int> channelNum;
std::vector<double> pulseH;
std::vector<double> digiTime;
std::vector<double> PCtime;
*/

/*
std::vector<std::vector<double> > pulseH(4, std::vector<double>(NTRIGS));
std::vector<std::vector<double> > digiTime(4, std::vector<double>(NTRIGS));
std::vector<std::vector<double> > PCtime(4, std::vector<double>(NTRIGS));
*/

std::vector<double> pulseH[8] = {std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>()};

std::vector<double> digiTime[8] = {std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>(),
				   std::vector<double>()};

std::vector<double> PCtime[8] = {std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>(),
				 std::vector<double>()};

std::vector<double> runStartTime[8] = {std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>(),
				       std::vector<double>()};


std::vector<int> rollover[8] = {std::vector<int>(),
				std::vector<int>(),
				std::vector<int>(),
				std::vector<int>(),
				std::vector<int>(),
				std::vector<int>(),
				std::vector<int>(),
				std::vector<int>()};


/****************************************************************************************************************/
/*														*/
/*					Digitizer Initialization Code						*/
/*														*/
/****************************************************************************************************************/

#include <CAENDigitizer.h>

 
#define MAXNB   1
#define MaxNChannels 8  
#define MAXNBITS 14
#include "Functions.h"

std::string path = "./";

bool params1=false;
bool params2=false;

bool first1=true;
bool first2=true;

int DCOffset = 0x1999;

void ClearHistos(void);
void QuitProgram(void);

/* The following variable is the type returned from most of CAENDigitizer
   library functions and is used to check if there was an error in function
   execution. For example:
     ret = CAEN_DGTZ_some_function(some_args);
     if(ret) printf("Some error"); */ 
CAEN_DGTZ_ErrorCode ret;
int iret;

/* Buffers to store the data. The memory must be allocated using the appropriate
   CAENDigitizer API functions (see below), so they must not be initialized here
   NB: you must use the right type for different DPP analysis (in this case PHA) */ 
char *buffer = NULL;					// readout buffer
CAEN_DGTZ_DPP_PHA_Event_t *Events[MaxNChannels];	// events buffer
CAEN_DGTZ_DPP_PHA_Waveforms_t *Waveform = NULL;		// waveforms buffer
  
/* The following variables will store the digitizer configuration parameters */ 
CAEN_DGTZ_DPP_PHA_Params_t DPPParams[MAXNB];
DigitizerParams_t Params[MAXNB];

/* Arrays for data analysis */ 
uint64_t PrevTime[MAXNB][MaxNChannels];
uint64_t ExtendedTT[MAXNB][MaxNChannels];
uint32_t * EHisto[MAXNB][MaxNChannels];		// Energy Histograms 
int ECnt[MAXNB][MaxNChannels];			// Total trigger count - doesn't get cleared each cycle
int PUCnt[MAXNB][MaxNChannels];			// Total pileup count - doesn't get cleared each cycle
int TrgCnt[MAXNB][MaxNChannels];		// Trigger count - gets cleared each cycle
int PurCnt[MAXNB][MaxNChannels];		// Pileup count - gets cleared each cycle

double Trg1Sec[MaxNChannels] = {0};
int IntegratedCount[MaxNChannels] = {0};
  
/* The following variable will be used to get an handler for the digitizer. The
   handler will be used for most of CAENDigitizer functions to identify the board */ 
int handle[MAXNB];

int USBhandle;

/* Other variables */ 
int i, b, ch, ev;
int Quit = 0;
int AcqRun = 0;
uint32_t AllocatedSize, BufferSize;
int Nb = 0;
int DoSaveWave[MAXNB][MaxNChannels];
int MajorNumber;
int BitMask = 0;
uint64_t CurrentTime, PrevRateTime, ElapsedTime, RunEndTime;
double RunStartTime;
uint32_t NumEvents[MaxNChannels];
CAEN_DGTZ_BoardInfo_t BoardInfo;

float enorm[MaxNChannels];

//FILE *out;	// file pointer to output file


extern "C" {
  long get_time ();
  void PrintInterface ();
  void Sleep (int t);
  int SaveHistogram (char *basename, int b, int ch, uint32_t * EHisto);
  int SaveDigitalProbe (int b, int ch, int trace, int size, uint8_t * WaveData);
  int SaveWaveform (int b, int ch, int trace, int size, int16_t * WaveData);
}



/****************************************************************************************************************/
/*                                                                                                              */
/*                                      Some helper methods                                                     */
/*                                                                                                              */
/****************************************************************************************************************/


//prints DPPParameters to file or console
void PrintParameters(FILE *metaData){

  int b=0;
  int ch=0;

  fprintf(metaData,"Decimation of Input Signal \t\t%d\n", DPPParams[b].decimation[ch]);
  fprintf(metaData,"Trigger Threshold \t\t\t%d\n", DPPParams[b].thr[ch]);
  fprintf(metaData,"Decay Time Constant (ns) \t\t%d\n", DPPParams[b].M[ch]);
  fprintf(metaData,"Trapezoid Rise Time (ns) \t\t%d\n", DPPParams[b].k[ch]);
  fprintf(metaData,"Trapezoid Flat Top (ns) \t\t%d\n", DPPParams[b].m[ch]);
  fprintf(metaData,"Flat top delay (peaking time) (ns) \t%d\n", DPPParams[b].ftd[ch]);
  fprintf(metaData,"Trigger Filter smoothing factor \t%d\n", DPPParams[b].a[ch]);
  fprintf(metaData,"Input Signal Rise time (ns) \t\t%d\n", DPPParams[b].b[ch]);
  fprintf(metaData,"Base Line Hold Off (ns) \t\t%d\n", DPPParams[b].blho[ch]);
  fprintf(metaData,"Trigger Hold Off (ns) \t\t\t%d\n", DPPParams[b].trgho[ch]);
  fprintf(metaData,"Peak Hold Off (ns) \t\t\t%d\n", DPPParams[b].pkho[ch]);
  fprintf(metaData,"Number of Samples for Baseline Mean \t%d\n", DPPParams[b].nsbl[ch]);
  fprintf(metaData,"Number of Samples for Peak Mean \t%d\n", DPPParams[b].nspk[ch]);
  fprintf(metaData,"Digital Probe Gain \t\t\t%d\n", DPPParams[b].dgain[ch]);
  fprintf(metaData,"Energy Normalization Factor \t\t%.2f\n", DPPParams[b].enf[ch]);
  fprintf(metaData,"Enable overlap rejection \t\t%d\n", DPPParams[b].otrej[ch]);
  fprintf(metaData,"ENABLE_RT_DISCR \t\t\t%d\n", DPPParams[b].trgwin[ch]);
  fprintf(metaData,"RT_DISCR_WINDOW \t\t\t%d\n", DPPParams[b].twwdt[ch]);
  
  fprintf(metaData,"\nDC offset\t\t\t\t%#05x\n", DCOffset);

  
}

//returns the time of the call.
double markTime() {
  gettimeofday(&timeMark,NULL);
  return (double)timeMark.tv_sec + (double)timeMark.tv_usec/1000000.; 
}


//writes the data to file
void writeNtuple(int k) {


  printf("Writing %d events to", m_numTrigs[k]);

  double fileStartName = runStartTime[k][0] + digiTime[k][0]*1e-8;

  //create a filename with the time of the first trigger in the filename
  std::string filename = path + Form("%s_%d_%F.root", "He3tubeData", k, fileStartName);

  char *y = new char[filename.length() + 1];
  std::strcpy(y, filename.c_str());
  
  printf(" %s\n", y);

  delete[] y;

  TString a = filename;
  TFile *outfile = new TFile(a,"RECREATE");

  int channel;
  double pHeight;  
  double dTime;
  double pcTime;
  //double runStart;
  int roll;

  //create a ttree and set the branches
  TTree *outtree = new TTree("tout", "tout");
  outtree->Branch("Channel", &channel, "Channel/I");
  outtree->Branch("DigiTime", &dTime, "DigiTime/D");          
  outtree->Branch("pulseHeight", &pHeight, "pulseHeight/D"); 
  outtree->Branch("PCtime", &pcTime, "PCtime/D");
  outtree->Branch("ExtendedTT", &roll, "ExtendedTT/I");

  for(int i=0; i<pulseH[k].size(); i++){
    channel = k;
    pHeight = pulseH[k][i];
    dTime = runStartTime[k][i]+digiTime[k][i]*1e-8;
    pcTime = PCtime[k][i];
    //runStart = runStartTime[k][i];
    roll = rollover[k][i];
    outtree->Fill();
  }

  //clear the vectors that hold the data
  printf("clearing vectors\n");
  pulseH[k].clear();
  digiTime[k].clear();
  PCtime[k].clear();
  runStartTime[k].clear();
  rollover[k].clear();
  
  //write to file.
  outtree->Write();
  outfile->Close();
}

FILE *fp;
std::stringstream ss;

void readIntegratedCountFile(){
  
  std::string fileBase = "he3Count_";
  int error=0;

  for(int i=0; i<8; i++){
    ss<<i;
    std::string filename = path+fileBase+ss.str()+".dat";
    ss.str(std::string());
    char *y = new char[filename.length() + 1];
    std::strcpy(y, filename.c_str());

    fp=fopen(y, "r");
    if (fp == NULL) 
      { 
	printf("error: could not open file %s\n", y);
	error=1;
      } else {
      int val;
      fscanf(fp, "%d", &val);
      printf("initial value from file: %d\n", val);
      fclose(fp);
      IntegratedCount[i] = val;
    }    
    delete[] y;
    if(error==1) IntegratedCount[i] = 0;
  }
}

void writeIntegratedCountFile(){

  std::string fileBase = "he3Count_";
  
  //printf("writing integrated counts to file\n");
  
  for(int i=0; i<8; i++){
    ss<<i;
    std::string filename = path+fileBase+ss.str()+".dat";
    ss.str(std::string());
    char *y = new char[filename.length() + 1];
    
    std::strcpy(y, filename.c_str());

    fp = fopen(y, "w");
    if (fp == NULL)
      {
      printf("error: could not open file %s\n", y);
      } else {
      fprintf(fp, "%d %u", IntegratedCount[i], (unsigned)time(NULL));
      fclose(fp);
    }
    delete[] y;
  }


}





int ProgramDigitizer (int handle, DigitizerParams_t Params, CAEN_DGTZ_DPP_PHA_Params_t DPPParams) 
{
  /* This function uses the CAENDigitizer API functions to perform the digitizer's initial configuration */ 
  int i, ret = 0;
    /* Reset the digitizer */ 
  ret |= CAEN_DGTZ_Reset (handle);
  if (ret)
  {
    printf ("ERROR: can't reset the digitizer.\n");
    return -1;
  }
  ret |= CAEN_DGTZ_WriteRegister (handle, 0x8000, 0x01000114);	// Channel Control Reg

  /* Set the DPP acquisition mode
     This setting affects the modes Mixed and List    (see CAEN_DGTZ_DPP_AcqMode_t definition for details)
     CAEN_DGTZ_DPP_SAVE_PARAM_EnergyOnly              Only energy is returned
     CAEN_DGTZ_DPP_SAVE_PARAM_TimeOnly                Only time is returned
     CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime           Both energy and time are returned
     CAEN_DGTZ_DPP_SAVE_PARAM_None                    No histogram data is returned */ 
  ret |= CAEN_DGTZ_SetDPPAcquisitionMode (handle, Params.AcqMode, CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime);
  
  // Set the digitizer acquisition mode (CAEN_DGTZ_SW_CONTROLLED or CAEN_DGTZ_S_IN_CONTROLLED)
  ret |= CAEN_DGTZ_SetAcquisitionMode (handle, CAEN_DGTZ_SW_CONTROLLED);
 
  // Set the number of samples for each waveform
  ret |= CAEN_DGTZ_SetRecordLength (handle, Params.RecordLength);

  // Set the I/O level (CAEN_DGTZ_IOLevel_NIM or CAEN_DGTZ_IOLevel_TTL)
  ret |= CAEN_DGTZ_SetIOLevel (handle, Params.IOlev);
  
  

  /* Set the digitizer's behaviour when an external trigger arrives:
     CAEN_DGTZ_TRGMODE_DISABLED: do nothing
     CAEN_DGTZ_TRGMODE_EXTOUT_ONLY: generate the Trigger Output signal
     CAEN_DGTZ_TRGMODE_ACQ_ONLY = generate acquisition trigger
     CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT = generate both Trigger Output and acquisition trigger
     see CAENDigitizer user manual, chapter "Trigger configuration" for details */ 
  ret |= CAEN_DGTZ_SetExtTriggerInputMode (handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);

  // Set the enabled channels
  ret |= CAEN_DGTZ_SetChannelEnableMask (handle, Params.ChannelMask);
  // Set how many events to accumulate in the board memory before being available for readout
  ret |= CAEN_DGTZ_SetDPPEventAggregation (handle, Params.EventAggr, 0);
  // Set the mode used to syncronize the acquisition between different boards.
  ret |= CAEN_DGTZ_SetRunSynchronizationMode (handle, CAEN_DGTZ_RUN_SYNC_Disabled);
  // Set the DPP specific parameters for the channels in the given channelMask
  ret |= CAEN_DGTZ_SetDPPParameters (handle, Params.ChannelMask, &DPPParams);

  for (i = 0; i < MaxNChannels; i++)
  {
    if (Params.ChannelMask & (1 << i))
    {  
      // Set a DC offset to the input signal to adapt it to digitizer's dynamic range
      //ret |= CAEN_DGTZ_SetChannelDCOffset (handle, i, 0x2666);
     //ret |= CAEN_DGTZ_SetChannelDCOffset (handle, i, DCOffset);
     //ret |= CAEN_DGTZ_SetChannelDCOffset (handle, i, 0x0ccc);
     ret |= CAEN_DGTZ_SetChannelDCOffset (handle, i, DCOffset);

      // Set the Pre-Trigger size (in samples)
      ret |= CAEN_DGTZ_SetDPPPreTriggerSize (handle, i, 200);
  
      // Set the polarity for the given channel (CAEN_DGTZ_PulsePolarityPositive or CAEN_DGTZ_PulsePolarityNegative)
      ret |= CAEN_DGTZ_SetChannelPulsePolarity (handle, i, Params.PulsePolarity);
    }
  }
  
  ret |= CAEN_DGTZ_SetDPP_PHA_VirtualProbe (handle,
					    CAEN_DGTZ_DPP_VIRTUALPROBE_DUAL,
				            CAEN_DGTZ_DPP_PHA_VIRTUALPROBE1_trapezoid,
				            CAEN_DGTZ_DPP_PHA_VIRTUALPROBE2_Input,
				            CAEN_DGTZ_DPP_PHA_DIGITAL_PROBE_Peaking);
  if (ret)
  {
    printf ("Warning: errors found during the programming of the digitizer.\nSome settings may not be executed\n");
    return ret;
  }
  else
  {
    return 0;
  }
}

bool DigiInit=false;
int InitializeDigitizer() 
{
  memset (DoSaveWave, 0, MAXNB * MaxNChannels * sizeof (int));

  for (i = 0; i < MAXNBITS; i++)
    BitMask |= 1 << i;		/* Create a bit mask based on number of bits of the board */
 

  /* *************************************************************************************** */ 
  /* Set Parameters                                                                          */ 
  /* *************************************************************************************** */ 
  memset (&Params, 0, MAXNB * sizeof (DigitizerParams_t));
  //memset (&DPPParams, 0, MAXNB * sizeof (CAEN_DGTZ_DPP_PHA_Params_t));
  for (b = 0; b < MAXNB; b++)
  {
    for (ch = 0; ch < MaxNChannels; ch++)
      EHisto[b][ch] = NULL;	//set all histograms pointers to NULL (we will allocate them later)
      
    /****************************\
    * Communication Parameters   *
    \****************************/ 
    // USB connection to V1718 bridge and access to the board with VME bus
    Params[b].LinkType = CAEN_DGTZ_USB;		// Link Type
    Params[b].VMEBaseAddress = 0x32100000;	// VME Base Address
    Params[b].IOlev = CAEN_DGTZ_IOLevel_TTL;
      
    /****************************\
    *  Acquisition parameters    *
    \****************************/ 
    //Params[b].AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope;	// CAEN_DGTZ_DPP_ACQ_MODE_List or CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope
    Params[b].AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_List;		// CAEN_DGTZ_DPP_ACQ_MODE_List or CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope
    Params[b].RecordLength = 1024;				// Num of samples of the waveforms (only for Oscilloscope mode)
    Params[b].ChannelMask = 0xFF;				// Channel enable mask - channels 0 -> 7
    //Params[b].ChannelMask = 0x1;				// Channel enable mask - channel 0 only               0x1 -> 1
    //Params[b].ChannelMask = 0x3;				// Channel enable mask - channels 0 and 1 only        0x3 -> 11
    //Params[b].ChannelMask = 0x7;				// Channel enable mask - channels 0, 1, and 2 only    0x7 -> 111
    //Params[b].ChannelMask = 0xF;                                // Channel enable mask - channels 0-3?                0xF -> 1111
    Params[b].EventAggr = 0;					// number of events in one aggregate (0=automatic)
    Params[b].PulsePolarity = CAEN_DGTZ_PulsePolarityNegative;	// Pulse Polarity (this parameter can be individual)
      
    /****************************\
    *      DPP parameters        *
    \****************************/ 
    
    if(!params1||!params2){
      printf("DPP parameters not set!\n");
      return 0;
	
    }
    
 
    // The energy (trapezoid height) is renormalized to make it (approximately...) independent
    // of the values chosen for k and M.
    for(i=0; i<MaxNChannels; i++) {
      if (Params[b].ChannelMask & (1<<i)) {
        int shf;
        float enf = DPPParams[b].enf[i];
        int k = DPPParams[b].k[i]/10;
        int M = DPPParams[b].M[i]/10;
        float ediv[MaxNChannels];
        for (shf=0; shf<64; shf++) {
          ediv[i] = (k*M)/(enf*(uint64_t)(1<<shf));
          if (ediv[i] < 2)
            break;
        }
        enorm[i] = ediv[i];
      }
    }  
  }
  /* *************************************************************************************** */ 
  /* Open the digitizer and read board information                                           */ 
  /* *************************************************************************************** */ 
  /* The following function is used to open the digitizer with the given connection parameters
     and get the handler to it */ 
  for (b = 0; b < MAXNB; b++)
  {
    /* IMPORTANT: The following function identifies the different boards with a system which may change
       for different connection methods (USB, Conet, ecc). Refer to CAENDigitizer user manual for more info.
       The following is for b boards connected to A2818 (or A3818) via opticalLink (or USB with A1718)
       in this case the boards are accessed throught VME bus, and you must specify the VME address of each board:
       Params[b].LinkType = CAEN_DGTZ_PCI_OpticalLink (CAEN_DGTZ_PCIE_OpticalLink for A3818 or CAEN_DGTZ_USB for A1718)
       Params[0].VMEBaseAddress = <0xXXXXXXXX> (address of first board) 
       Params[1].VMEBaseAddress = <0xYYYYYYYY> (address of second board) 
       etc */ 
    ret = CAEN_DGTZ_OpenDigitizer (Params[b].LinkType, 0, 0, Params[b].VMEBaseAddress, &handle[b]);
    if (ret)
    {
      printf ("Can't open digitizer\n");
      return 0;
      //QuitProgram();
    }

    /* Once we have the handler to the digitizer, we use it to call the other functions */ 
    ret = CAEN_DGTZ_GetInfo (handle[b], &BoardInfo);
    if (ret)
    {
      printf ("Can't read board info\n");
      QuitProgram();
    }
    printf ("\nConnected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, b);
    printf ("ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
    printf ("AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);

    //DigiInit=true;
      
    /* Check firmware revision (only DPP firmwares can be used with this Demo) */ 
    sscanf (BoardInfo.AMC_FirmwareRel, "%d", &MajorNumber);
    if (MajorNumber != 128)
    {
      printf ("This digitizer has not a DPP-PHA firmware\n");
      QuitProgram();
    }
  }
  
  /* *************************************************************************************** */ 
  /* Program the digitizer (see function ProgramDigitizer)                                   */ 
  /* *************************************************************************************** */ 
  for (b = 0; b < MAXNB; b++)
  {
    iret = ProgramDigitizer (handle[b], Params[b], DPPParams[b]);
    if (iret)
    {
      printf ("Failed to program the digitizer\n");
      QuitProgram();
    }
  }
  
  /* WARNING: The mallocs MUST be done after the digitizer programming,
     because the following functions needs to know the digitizer configuration
     to allocate the right memory amount */

  /* Allocate memory for the readout buffer */ 
  ret = CAEN_DGTZ_MallocReadoutBuffer (handle[0], &buffer, &AllocatedSize);
  
  /* Allocate memory for the events */ 
  iret |= CAEN_DGTZ_MallocDPPEvents (handle[0], (void**) Events, &AllocatedSize);
  
  /* Allocate memory for the waveforms */ 
  iret |= CAEN_DGTZ_MallocDPPWaveforms (handle[0], (void**) &Waveform, &AllocatedSize);
  if (iret)
  {
    printf ("Can't allocate memory buffers\n");
    QuitProgram();
  }

  //out = fopen("data.dat", "w");	// open output file

  DigiInit=true;
  return 0;
}


/****************************************************************************************************************/
/*														*/
/*						Readout Loop							*/
/*														*/
/****************************************************************************************************************/
void readoutData (void *ptr) {

  printf ("Starting readout loop...\n");

  
  double oldReadTime=0;
  double newReadTime=0;

  ClearHistos();

  for (b = 0; b < MAXNB; b++)
  {
    // Start Acquisition
    // NB: the acquisition for each board starts when the following line is executed
    // so in general the acquisition does NOT starts syncronously for different boards
    CAEN_DGTZ_SWStartAcquisition (handle[b]);
    RunStartTime = markTime();

    printf ("Acquisition Started for Board %d\n", b);
  }
  //  RunStartTime = markTime();

  while (true) {
    while (mq1->pending()==0) {

      newReadTime = markTime();

      if (! AcqRun) {
        sleep (2);
        continue;
      }

      uint32_t firsttime, lasttime;
      int lpnum;

      /* Read data from the boards */ 
      for (b = 0; b < MAXNB; b++)
      {  
        /* Read data from the board */ 
        ret = CAEN_DGTZ_ReadData (handle[b], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer, &BufferSize);
	//ret = CAEN_DGTZ_ReadData (handle[b], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_2eSST, buffer, &BufferSize);
	

	oldReadTime = newReadTime; 
	newReadTime = markTime();


	if (ret)
        {
          printf ("Readout Error\n");
          QuitProgram();
        }
        if (BufferSize == 0) {	// no triggers...
          continue;
        }
        Nb += BufferSize;

        //ret = DataConsistencyCheck((uint32_t *)buffer, BufferSize/4);
        iret |= CAEN_DGTZ_GetDPPEvents (handle[b], buffer, BufferSize, (void**) Events, NumEvents);
	if (ret)
        {
          printf ("Data Error: %d\n", ret);
          QuitProgram();
        }


        /* Analyze data */ 
        for (ch = 0; ch < MaxNChannels; ch++)
        {
          if (!(Params[b].ChannelMask & (1 << ch))) continue;

	  firsttime=1;
	  lasttime=0;
	  lpnum=0;

          /* Update Histograms */ 
          for (ev = 0; ev < NumEvents[ch]; ev++)
          {
	    TrgCnt[b][ch]++;  
	    //Trg1Sec[ch]++;
	    
	    //printf("channel  number of triggers: %d\n", Trg1Sec[ch]);

	    m_numTrigs[ch]++;
	    IntegratedCount[ch]++;
            /* Time Tag */ 
            /* Events[ch][ev].TimeTag is in units of 10 ns, and starts at 0 at the beginning of the run. */
            /* At some point, it must eventually roll over... Maybe ExtendedTT takes care of this??? */
            /* But if no triggers for a VERY long time, could ExtendedTT miss an entire TimeTag cycle??? */
            /* Not to worry - TimeTag is a uint64, so it won't roll over for 5845 years!!! */
            /* (assuming the FPGA *also* allocates 64 bits...) so ExtendedTT is probably useless. */
            if (Events[ch][ev].TimeTag < PrevTime[b][ch]) {	// TimeTag rolled over...
              ExtendedTT[b][ch]++;
            }
            PrevTime[b][ch] = Events[ch][ev].TimeTag;
	    
	    if(lpnum==0) firsttime=Events[ch][ev].TimeTag;
	    else if(lpnum==NumEvents[ch]-1) lasttime=Events[ch][ev].TimeTag;

            /* Energy */ 
            if (Events[ch][ev].Energy > 0)	// Fill the histograms
            {
              uint32_t myEnergy = (uint32_t) ((Events[ch][ev].Energy) & ((1<<(MAXNBITS+1))-1));
              myEnergy = myEnergy / enorm[ch];
	      ECnt[b][ch]++;

	      //channelNum.push_back(ch);
	      pulseH[ch].push_back((double)myEnergy);

	      //the digitizer clock starts at RunStartTime, TimeTag is in units of 10ns since the run began
	      digiTime[ch].push_back((double)Events[ch][ev].TimeTag);
	      PCtime[ch].push_back((double)oldReadTime);
	      rollover[ch].push_back(ExtendedTT[b][ch]);
	      runStartTime[ch].push_back(RunStartTime);
	    }
            else	// PileUp
            {
	      
	      //channelNum.push_back(ch);
              pulseH[ch].push_back(-10);   //energy value of -10 means a pileup has occured.
              digiTime[ch].push_back((double)Events[ch][ev].TimeTag);
	      PCtime[ch].push_back((double)oldReadTime);
	      runStartTime[ch].push_back(RunStartTime);
	      rollover[ch].push_back(ExtendedTT[b][ch]);


              PUCnt[b][ch]++;
              PurCnt[b][ch]++;

	    }
	    
	    //when NTRIGS events have happened, save the data to file
	    if(m_numTrigs[ch]==NTRIGS) {               
	      writeNtuple(ch);                          // Write Ntuples
	      m_numTrigs[ch]=0;                           // reset the trigger count
	    }
	    lpnum++;
          }	// loop on events
	  if(NumEvents[ch]!=0) Trg1Sec[ch] = (double)NumEvents[ch]/((double)(lasttime-firsttime)*1e-8);
	  IntegratedCount[ch] = IntegratedCount[ch]-numSwTrigs;
	  numSwTrigs=0;
	}	// loop on channels
      }		// loop on boards
    }		// end of readout loop

    int receive=0;
    mq1->receive(&receive, sizeof(int));	// get the message from the queue
    printf ("received %d\n", receive);		// print out message

    switch (receive) {

      case 2: {					// clear histograms
        printf("Clearing histograms...\n");
        ClearHistos();
      } break;

      case 3: {					// update histograms
        printf("Updating histograms...\n");
        /*
	  for (b = 0; b < MAXNB; b++)
          for (ch = 0; ch < MaxNChannels; ch++)
            if (ECnt[b][ch] != 0)              
	    SaveHistogram ("Histo", b, ch, EHisto[b][ch]);	// Save Histograms to file for each board
	*/
      } break;

      case 4: {					// pause daq
        printf ("pausing daq...\n");
        AcqRun = 0;
      } break;

      case 5: {
        printf ("resuming daq...\n");
        AcqRun = 1;
      } break;

      case 6: {					// reset myThread and exit readout loop
        printf("Acquisition Ended\n");
	for(int k=0; k<4; k++){
	  if(m_numTrigs[k]!=0) writeNtuple(k);  //write data to file.
	  
	}
	writeIntegratedCountFile();
	
        QuitProgram();
        plock->lock();
        myThread=0; 
        plock->unlock();
        AcqRun = 0;
	isRunning=false;
        return;
      } break;

      default: {
        printf ("undefined signal...\n");
      } break;

    }		// end of switch (recieve)
  }		// end of while (true)


}

/****************************************************************************************************************/
/*														*/
/*					Digitizer Utility Routines						*/
/*														*/
/****************************************************************************************************************/

void QuitProgram () 
{
  // stop the acquisition, close the device and free the buffers
  for (b = 0; b < MAXNB; b++)
  {
    CAEN_DGTZ_SWStopAcquisition (handle[b]);    
    CAEN_DGTZ_CloseDigitizer (handle[b]);
    for (ch = 0; ch < MaxNChannels; ch++)
      if (EHisto[b][ch] != NULL)
        free (EHisto[b][ch]);
  }
  CAEN_DGTZ_FreeReadoutBuffer (&buffer);
  CAEN_DGTZ_FreeDPPEvents (handle[0], (void**) Events);
  CAEN_DGTZ_FreeDPPWaveforms (handle[0], Waveform);

  RunEndTime = get_time ();
  return;
}

void ClearHistos ()
{
  // Clear Histograms and counters
  for (b = 0; b < MAXNB; b++)
  {
    for (ch = 0; ch < MaxNChannels; ch++)
    {
      EHisto[b][ch] = (uint32_t *) malloc ((1 << MAXNBITS) * sizeof (uint32_t));
      memset (EHisto[b][ch], 0, (1 << MAXNBITS) * sizeof (uint32_t));
      TrgCnt[b][ch] = 0;
      ECnt[b][ch] = 0;
      PUCnt[b][ch] = 0;
      PrevTime[b][ch] = 0;
      ExtendedTT[b][ch] = 0;
      PurCnt[b][ch] = 0;
    }
  }
}

/****************************************************************************************************************/
/*				  	        EPICS code							*/
/*														*/
/****************************************************************************************************************/

//uses the subRecord to send commands to the digitizer, such as initalize, start, stop, etc...
//setting of digitizer parameters will likely be implemented here
long commandRoutine(struct subRecord *psub) {
  switch ((int)psub->a) {
    //could be moved to main method?
  case 0: {
    message=6;
      plock->lock();
      if(isRunning) {
	mq1->send((void *) &message, sizeof(int));
	isRunning=false;
	psub->l=0;
      } else {
	printf("Aquisition not started yet!");
	//psub->a=6;
      }
      plock->unlock();
    } break;					

  case 1: {						// start readout loop                  
      printf ("Initializing the digitizer...\n");
      InitializeDigitizer();
      if(!DigiInit){
	psub->l=0;
	printf("Digitizer not initalized properly\n");
	break;	
      }
      
      plock->lock();// lock resource
      if (myThread==0){					// no thread in existance
        AcqRun = 1;
	isRunning=true;
        myThread=epicsThreadCreate("myThread", epicsThreadPriorityMedium, epicsThreadStackSmall, readoutData, NULL);
	psub->l=1;
      }
      else printf( "There's already a thread!\n");	//if there is a thread, do nothing and print message
      plock->unlock();					//release resource
    } break;

    case 2: {						// clear histograms
      message=2;
      mq1->send ((void *) &message, sizeof(int));	//sends message to epicsMessageQueue (mq1)
    } break;

    case 3: {						// update histograms
      message=3;
      mq1->send ((void *) &message, sizeof(int)); 
    } break;

    case 4: {						// pause daq
      message=4;
      mq1->send((void *) &message, sizeof(int));
    } break;

    case 5: {						// resume daq
      message=5;
      mq1->send((void *) &message, sizeof(int));
    } break;

    case 6: {
    } break;
  
    default: {
      printf ("undefined signal...\n");
    } break;
  }

  printf("Digitizer control record called\n");

  return(0);
};

static long paramRoutine1(struct subRecord *psub) {

  if(!params1&&!params2) memset (&DPPParams, 0, MAXNB * sizeof (CAEN_DGTZ_DPP_PHA_Params_t));
  
  if(psub->a!=-1){
    for(int b=0; b<MAXNB; b++){
      for (ch = 0; ch < MaxNChannels; ch++){
	DPPParams[b].decimation[ch] = psub->b;	// Decimation of Input Signal
	DPPParams[b].thr[ch] = psub->c;		// Trigger Threshold
	DPPParams[b].M[ch] = psub->d;		// Decay Time Constant (ns)
	DPPParams[b].k[ch] = psub->e;		// Trapezoid Rise Time (ns)
	DPPParams[b].m[ch] = psub->f;		// Trapezoid Flat Top (ns)
	DPPParams[b].ftd[ch] = psub->g;		// Flat top delay (peaking time) (ns)
	DPPParams[b].a[ch] = psub->h;		// Trigger Filter smoothing factor
	DPPParams[b].b[ch] = psub->i;		// Input Signal Rise time (ns)
	DPPParams[b].blho[ch] = psub->j;	// Base Line Hold Off (ns)
	DPPParams[b].trgho[ch] = psub->k;	// Trigger Hold Off (ns)
	DCOffset = psub->l;
      }
    }
    params1=true;
  }

  //print out values of parameters
  if(!first1){
    PrintParameters(stdout);
  } else{
    printf("First set of DPP parameters set\n");
    psub->a=-1;
  }
  first1=false;

  return(0);
}

static long paramRoutine2(struct subRecord *psub) {

  if(!params1&&!params2) memset (&DPPParams, 0, MAXNB * sizeof (CAEN_DGTZ_DPP_PHA_Params_t));

  if(psub->a!=-1){
    for(int b=0; b<MAXNB; b++){
      for (ch = 0; ch < MaxNChannels; ch++){

	DPPParams[b].nsbl[ch] = psub->b;		// Number of Samples for Baseline Mean; 3 = bx10 = 64 samples
	DPPParams[b].nspk[ch] = psub->c;		// Number of Samples for Peak Mean
	DPPParams[b].dgain[ch] = psub->d;		// Digital Probe Gain
	DPPParams[b].enf[ch] = psub->e;	        	// Energy Normalization Factor
	DPPParams[b].otrej[ch] = psub->f;		// Enable overlap rejection
	DPPParams[b].trgwin[ch] = psub->g;		// ENABLE_RT_DISCR
	DPPParams[b].twwdt[ch] = psub->h;		// RT_DISCR_WINDOW
	DPPParams[b].pkho[ch] = psub->i;                // Peak hold off (ns)

      }
    }
    params2=true;
  }

  NTRIGS = (int)psub->j;

  PrintParameters(stdout);
  
  printf("Saving ntuple after %d events\n", NTRIGS);
  
  if(first2) psub->a=-1;
  first2=false;

  return(0);

}

static long setPath(struct subRecord *psub) {
  path=psub->desc;
  std::string message = "Data path set to "+path+"\n";
  printf(message.c_str());

  readIntegratedCountFile();
  
  return 0;
}

static long getRates(struct subRecord *psub) {

  psub->a=Trg1Sec[0];
  psub->b=Trg1Sec[1];
  psub->c=Trg1Sec[2];
  psub->d=Trg1Sec[3];
  psub->e=Trg1Sec[4];
  psub->f=Trg1Sec[5];
  psub->g=Trg1Sec[6];
  psub->h=Trg1Sec[7];

  //pileups
  psub->e=PUCnt[0][0];  
  psub->f=PUCnt[0][1];
  psub->g=PUCnt[0][2];
  psub->h=PUCnt[0][3];

  if(isRunning){
    psub->i=1;
  }else {
    psub->i=0;
  }

  return 0;

}

static long getPileup(struct subRecord *psub) {

  //pileups
  psub->a=PUCnt[0][0];  
  psub->b=PUCnt[0][1];
  psub->c=PUCnt[0][2];
  psub->d=PUCnt[0][3];
  psub->e=PUCnt[0][4];  
  psub->f=PUCnt[0][5];
  psub->g=PUCnt[0][6];
  psub->h=PUCnt[0][7];

  return 0;
}

static long getIntegratedHits(struct subRecord *psub) {

  psub->a = IntegratedCount[0];
  psub->b = IntegratedCount[1];
  psub->c = IntegratedCount[2];
  psub->d = IntegratedCount[3];

  if(isRunning) writeIntegratedCountFile();


  return 0;
}

static long sendSWtrigger(struct subRecord *psub){

  if(isRunning){
    if(psub->a==1) printf("software trigger sent at %f\n", markTime());
    for(int b=0; b<MAXNB; b++){
      CAEN_DGTZ_SendSWtrigger(handle[b]);
      numSwTrigs++;
    }
  }
  return 0;
}

static long saveData(struct subRecord *psub){
  printf("Saving buffered data due to EPICS request\n");
  if(isRunning){
    //pause DAQ
    message=4;
    mq1->send((void *) &message, sizeof(int));      
    while(AcqRun==0) sleep (1);
    
    //write data
    for(int k=0; k<8; k++){
      if(m_numTrigs[k]!=0) writeNtuple(k); 
      m_numTrigs[k]=0;
    }
    //resume DAQ
    message=5;
    mq1->send((void *) &message, sizeof(int));

  }
  return 0;
}


//register routines
epicsRegisterFunction(commandRoutine);  
epicsRegisterFunction(paramRoutine1);
epicsRegisterFunction(paramRoutine2);
epicsRegisterFunction(setPath);
epicsRegisterFunction(getRates);
epicsRegisterFunction(getPileup);
epicsRegisterFunction(getIntegratedHits);
epicsRegisterFunction(sendSWtrigger);
epicsRegisterFunction(saveData);
//note that these routines must also be registered in he3DAQApp/src/he3DAQcroutines.dbd



//main method, self explanatory
int main(int argc,char *argv[])
{

  mq1 = new epicsMessageQueue(10,20);  //initalize the message queue. 
  Data=new double[100];           //intitialize the array
  for(int i=0; i<100; i++){            //set all values to 0
    Data[i]=0;
  }
  

  plock = new epicsMutex;            //initalize mutex.
  if(plock==0){
    printf(" plock is null :(\n");
    return 1;
  }

  if(argc>=2) {    
    iocsh(argv[1]);
    epicsThreadSleep(.2);
  }
  
  iocsh(NULL);

  message=6;
  if(isRunning) {
    mq1->send((void *) &message, sizeof(int));
    isRunning=false;
  }
  
  epicsExit(0);

  return(0);


}
