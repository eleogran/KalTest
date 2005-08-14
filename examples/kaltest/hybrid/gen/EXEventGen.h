#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "TKalDetCradle.h"
#include "THelicalTrack.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt = 1.);
   void          Swim(THelicalTrack &heltrk);

private:
   TKalDetCradle *fCradlePtr;
   TObjArray     *fHitBufPtr;

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif