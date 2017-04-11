TOP=../..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

#=============================
# Build the IOC application

PROD_IOC = he3DAQ
# he3DAQ.dbd will be created and installed
DBD += he3DAQ.dbd

# he3DAQ.dbd will be made up from these files:
he3DAQ_DBD += base.dbd
he3DAQ_DBD += he3DAQcroutines.dbd

# Include dbd files from all support applications:
#he3DAQ_DBD += xxx.dbd

#add root libraries   					srdejong mods here

ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBDIR   := $(shell $(ROOTCONFIG) --libdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --noldflags --noauxlibs --libs | sed s/-l//g)

USR_CXXFLAGS += -Wall -O2 $(ROOTCFLAGS)

# Add all the support libraries needed by this IOC	pof mods here
SHRLIB_SEARCH_DIRS += /usr/lib
SHRLIB_SEARCH_DIRS += $(ROOTLIBDIR)

he3DAQ_LIBS += CAENDigitizer
he3DAQ_LIBS += $(ROOTLIBS)

# he3DAQ_registerRecordDeviceDriver.cpp derives from he3DAQ.dbd
he3DAQ_SRCS += he3DAQ_registerRecordDeviceDriver.cpp

# Build the main IOC entry point on workstation OSs.	pof mods here
he3DAQ_SRCS_DEFAULT += he3DAQMain.cpp
he3DAQ_SRCS_DEFAULT += Functions.c
he3DAQ_SRCS_DEFAULT += keyb.c
he3DAQ_SRCS_vxWorks += -nil-

# Finally link to the EPICS Base libraries
he3DAQ_LIBS += $(EPICS_BASE_IOC_LIBS)

#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

