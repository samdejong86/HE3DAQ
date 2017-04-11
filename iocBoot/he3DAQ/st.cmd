

#cd ${TOP}

## Register all support components
dbLoadDatabase "dbd/he3DAQ.dbd"
he3DAQ_registerRecordDeviceDriver pdbbase


## Load record instances
#ld < db/calcRecord.o
dbLoadRecords("db/he3DAQrecords.db","SubDet=BEAST:DAQ:HE3,user=samHost")

#cd ${TOP}/iocBoot/${IOC}

#ld < db/calcRecord.o
#dbLoadDatabase("cRecord.dbd")


iocInit

## Start any sequence programs
dbpf BEAST:DAQ:HE3:DIGI:PARAMS1.A 0
dbpf BEAST:DAQ:HE3:DIGI:PARAMS2.A 0
dbpf BEAST:DAQ:HE3:SETPATH.A 0
dbgrep "*"
#seq sncxxx,"user=samHost"