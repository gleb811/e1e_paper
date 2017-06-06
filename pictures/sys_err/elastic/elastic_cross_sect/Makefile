
CXX           = g++ -Wno-write-strings -Wno-pragmas
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c -L/apps/root/PRO/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lm -ldl -rdynamic
ROOTGLIBS    := -L/apps/root/PRO/root/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lm -ldl -rdynamic
ROOTINCLUDE  := -I$(shell root-config --incdir)



all:
	make Dictationarys
	make h10tot21

h10tot21: h10tot21.$(ObjSuf) MyMainFrame.$(ObjSuf) MyMainFrameDict.$(ObjSuf) macro.$(ObjSuf)  macroDict.$(ObjSuf) mom_corr.$(ObjSuf) mom_corrDict.$(ObjSuf) elaslib.o 
	$(CXX) -g -o $@ $^ $(ROOTGLIBS) 


	
Dictationarys:	
	rootcint -f MyMainFrameDict.cxx -c -I`root-config --incdir` MyMainFrame.h 
	rootcint -f macroDict.cxx -c -I`root-config --incdir` macro.h
	rootcint -f mom_corrDict.cxx -c -I`root-config --incdir` mom_corr.h
	g77 -c -fPIC elaslib.f

#	rootcint -f $(TARGETCINT) -c -I$(ROOTSYS)/include $(TARGET).h

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) -g -c $(ROOTINCLUDE) -c $<

clean:
	rm -f *.o
	rm -f *Dict.*
	rm -f G__*
	rm h10tot21
	

