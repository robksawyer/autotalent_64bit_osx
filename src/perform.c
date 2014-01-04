  double* pfInput;
  double* pfOutput;
  float fAmount;
  float fSmooth;
  int iNotes[12];
  int iPitch2Note[12];
  int iNote2Pitch[12];
  int numNotes;
  float fTune;
  float fFixed;
  float fPull;
  float fShift;
  int iScwarp;
  float fLfoamp;
  float fLforate;
  float fLfoshape;
  float fLfosymm;
  int iLfoquant;
  int iFcorr;
  float fFwarp;
  float fMix;
//  Autotalent* psAutotalent;
  unsigned long lSampleIndex;

  long int N;
  long int Nf;
  long int fs;
  float pmin;
  float pmax;
  unsigned long nmin;
  unsigned long nmax;

  long int ti;
  long int ti2;
  long int ti3;
  long int ti4;
  float tf;
  float tf2;

  // Variables for cubic spline interpolator
  float indd;
  int ind0;
  int ind1;
  int ind2;
  int ind3;
  float vald;
  float val0;
  float val1;
  float val2;
  float val3;

  int lowersnap;
  int uppersnap;
  float lfoval;

  float pperiod;
  float inpitch;
  float conf;
  float outpitch;
  float aref;
  float fa;
  float fb;
  float fc;
  float fk;
  float flamb;
  float frlamb;
  float falph;
  float foma;
  float f1resp;
  float f0resp;
  float flpa;
  int ford;
  //psAutotalent = (Autotalent *)Instance;

  pfInput = inputs[0];
  pfOutput = outputs[0];
  fAmount = (float) *(m_pfAmount);
  fSmooth = (float) *(m_pfSmooth) * 0.8; // Scales max to a more reasonable value
  fTune = (float) *(m_pfTune);
  iNotes[0] = (int) *(m_pfA);
  iNotes[1] = (int) *(m_pfBb);
  iNotes[2] = (int) *(m_pfB);
  iNotes[3] = (int) *(m_pfC);
  iNotes[4] = (int) *(m_pfDb);
  iNotes[5] = (int) *(m_pfD);
  iNotes[6] = (int) *(m_pfEb);
  iNotes[7] = (int) *(m_pfE);
  iNotes[8] = (int) *(m_pfF);
  iNotes[9] = (int) *(m_pfGb);
  iNotes[10] = (int) *(m_pfG);
  iNotes[11] = (int) *(m_pfAb);
  fFixed = (float) *(m_pfFixed);
  fPull = (float) *(m_pfPull);
  fShift = (float) *(m_pfShift);
  iScwarp = (int) *(m_pfScwarp);
  fLfoamp = (float) *(m_pfLfoamp);
  fLforate = (float) *(m_pfLforate);
  fLfoshape = (float) *(m_pfLfoshape);
  fLfosymm = (float) *(m_pfLfosymm);
  iLfoquant = (int) *(m_pfLfoquant);
  iFcorr = (int) *(m_pfFcorr);
  fFwarp = (float) *(m_pfFwarp);
  fMix = (float) *(m_pfMix);

  // Some logic for the semitone->scale and scale->semitone conversion
  // If no notes are selected as being in the scale, instead snap to all notes
  ti2 = 0;
  for (ti=0; ti<12; ti++) {
    if (iNotes[ti]>=0) {
      iPitch2Note[ti] = ti2;
      iNote2Pitch[ti2] = ti;
      ti2 = ti2 + 1;
    }
    else {
      iPitch2Note[ti] = -1;
    }
  }
  numNotes = ti2;
  while (ti2<12) {
    iNote2Pitch[ti2] = -1;
    ti2 = ti2 + 1;
  }
  if (numNotes==0) {
    for (ti=0; ti<12; ti++) {
      iNotes[ti] = 1;
      iPitch2Note[ti] = ti;
      iNote2Pitch[ti] = ti;
    }
    numNotes = 12;
  }
  iScwarp = (iScwarp + numNotes*5)%numNotes;

  ford = mford;
  falph = mfalph;
  foma = (float)1 - falph;
  flpa = mflpa;
  flamb = mflamb;
  tf = pow((float)2,fFwarp/2)*(1+flamb)/(1-flamb);
  frlamb = (tf - 1)/(tf + 1);

  aref = (float)fTune;

  N = cbsize;
  Nf = corrsize;
  fs = mfs;

  pmax = psAutotalent->pmax;
  pmin = psAutotalent->pmin;
  nmax = psAutotalent->nmax;
  nmin = psAutotalent->nmin;

  aref = psAutotalent->aref;
  pperiod = psAutotalent->pmax;
  inpitch = psAutotalent->inpitch;
  conf = psAutotalent->conf;
  outpitch = psAutotalent->outpitch;


  /*******************
   *  MAIN DSP LOOP  *
   *******************/
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++)  {
    
    // load data into circular buffer
    tf = (float) *(pfInput++);
    ti4 = psAutotalent->cbiwr;
    psAutotalent->cbi[ti4] = tf;

    if (iFcorr>=1) {
      // Somewhat experimental formant corrector
      //  formants are removed using an adaptive pre-filter and
      //  re-introduced after pitch manipulation using post-filter
      // tf is signal input
      fa = tf - psAutotalent->fhp; // highpass pre-emphasis filter
      psAutotalent->fhp = tf;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	psAutotalent->fsig[ti] = fa*fa*foma + psAutotalent->fsig[ti]*falph;
	fc = (fb-psAutotalent->fc[ti])*flamb + psAutotalent->fb[ti];
	psAutotalent->fc[ti] = fc;
	psAutotalent->fb[ti] = fb;
	fk = fa*fc*foma + psAutotalent->fk[ti]*falph;
	psAutotalent->fk[ti] = fk;
	tf = fk/(psAutotalent->fsig[ti] + 0.000001);
	tf = tf*foma + psAutotalent->fsmooth[ti]*falph;
	psAutotalent->fsmooth[ti] = tf;
	psAutotalent->fbuff[ti][ti4] = tf;
	fb = fc - tf*fa;
	fa = fa - tf*fc;
      }
      psAutotalent->cbf[ti4] = fa;
      // Now hopefully the formants are reduced
      // More formant correction code at the end of the DSP loop
    }
    else {
      psAutotalent->cbf[ti4] = tf;
    }


    // Input write pointer logic
    psAutotalent->cbiwr++;
    if (psAutotalent->cbiwr >= N) {
      psAutotalent->cbiwr = 0;
    }


    // ********************
    // * Low-rate section *
    // ********************

    // Every N/noverlap samples, run pitch estimation / manipulation code
    if ((psAutotalent->cbiwr)%(N/psAutotalent->noverlap) == 0) {

      // ---- Obtain autocovariance ----

      // Window and fill FFT buffer
      ti2 = psAutotalent->cbiwr;
      for (ti=0; ti<N; ti++) {
	psAutotalent->ffttime[ti] = (float)(psAutotalent->cbi[(ti2-ti+N)%N]*psAutotalent->cbwindow[ti]);
      }

      // Calculate FFT
      fft_forward(psAutotalent->fmembvars, psAutotalent->ffttime, psAutotalent->fftfreqre, psAutotalent->fftfreqim);

      // Remove DC
      psAutotalent->fftfreqre[0] = 0;
      psAutotalent->fftfreqim[0] = 0;

      // Take magnitude squared
      for (ti=1; ti<Nf; ti++) {
	psAutotalent->fftfreqre[ti] = (psAutotalent->fftfreqre[ti])*(psAutotalent->fftfreqre[ti]) + (psAutotalent->fftfreqim[ti])*(psAutotalent->fftfreqim[ti]);
	psAutotalent->fftfreqim[ti] = 0;
      }

      // Calculate IFFT
      fft_inverse(psAutotalent->fmembvars, psAutotalent->fftfreqre, psAutotalent->fftfreqim, psAutotalent->ffttime);

      // Normalize
      tf = (float)1/psAutotalent->ffttime[0];
      for (ti=1; ti<N; ti++) {
	psAutotalent->ffttime[ti] = psAutotalent->ffttime[ti] * tf;
      }
      psAutotalent->ffttime[0] = 1;

      //  ---- END Obtain autocovariance ----


      //  ---- Calculate pitch and confidence ----

      // Calculate pitch period
      //   Pitch period is determined by the location of the max (biased)
      //     peak within a given range
      //   Confidence is determined by the corresponding unbiased height
      tf2 = 0;
      pperiod = pmin;
      for (ti=nmin; ti<nmax; ti++) {
	ti2 = ti-1;
	ti3 = ti+1;
	if (ti2<0) {
	  ti2 = 0;
	}
	if (ti3>Nf) {
	  ti3 = Nf;
	}
	tf = psAutotalent->ffttime[ti];

	if (tf>psAutotalent->ffttime[ti2] && tf>=psAutotalent->ffttime[ti3] && tf>tf2) {
	  tf2 = tf;
	  ti4 = ti;
	}
      }
      if (tf2>0) {
	conf = tf2*psAutotalent->acwinv[ti4];
	if (ti4>0 && ti4<Nf) {
	  // Find the center of mass in the vicinity of the detected peak
	  tf = psAutotalent->ffttime[ti4-1]*(ti4-1);
	  tf = tf + psAutotalent->ffttime[ti4]*(ti4);
	  tf = tf + psAutotalent->ffttime[ti4+1]*(ti4+1);
	  tf = tf/(psAutotalent->ffttime[ti4-1] + psAutotalent->ffttime[ti4] + psAutotalent->ffttime[ti4+1]);
	  pperiod = tf/fs;
	}
	else {
	  pperiod = (float)ti4/fs;
	}
      }

      // Convert to semitones
      tf = (float) -12*log10((float)aref*pperiod)*L2SC;
      if (conf>=psAutotalent->vthresh) {
	inpitch = tf;
	psAutotalent->inpitch = tf; // update pitch only if voiced
      }
      psAutotalent->conf = conf;

      *(psAutotalent->m_pfPitch) = (LADSPA_Data) inpitch;
      *(psAutotalent->m_pfConf) = (LADSPA_Data) conf;

      //  ---- END Calculate pitch and confidence ----


      //  ---- Modify pitch in all kinds of ways! ----

      outpitch = inpitch;

      // Pull to fixed pitch
      outpitch = (1-fPull)*outpitch + fPull*fFixed;

      // -- Convert from semitones to scale notes --
      ti = (int)(outpitch/12 + 32) - 32; // octave
      tf = outpitch - ti*12; // semitone in octave
      ti2 = (int)tf;
      ti3 = ti2 + 1;
      // a little bit of pitch correction logic, since it's a convenient place for it
      if (iNotes[ti2%12]<0 || iNotes[ti3%12]<0) { // if between 2 notes that are more than a semitone apart
	lowersnap = 1;
	uppersnap = 1;
      }
      else {
	lowersnap = 0;
	uppersnap = 0;
	if (iNotes[ti2%12]==1) { // if specified by user
	  lowersnap = 1;
	}
	if (iNotes[ti3%12]==1) { // if specified by user
	  uppersnap = 1;
	}
      }
      // (back to the semitone->scale conversion)
      // finding next lower pitch in scale
      while (iNotes[(ti2+12)%12]<0) {
      	ti2 = ti2 - 1;
      }
      // finding next higher pitch in scale
      while (iNotes[ti3%12]<0) {
      	ti3 = ti3 + 1;
      }
      tf = (tf-ti2)/(ti3-ti2) + iPitch2Note[(ti2+12)%12];
      if (ti2<0) {
      	tf = tf - numNotes;
      }
      outpitch = tf + numNotes*ti;
      // -- Done converting to scale notes --

      // The actual pitch correction
      ti = (int)(outpitch+128) - 128;
      tf = outpitch - ti - 0.5;
      ti2 = ti3-ti2;
      if (ti2>2) { // if more than 2 semitones apart, put a 2-semitone-like transition halfway between
	tf2 = (float)ti2/2;
      }
      else {
	tf2 = (float)1;
      }
      if (fSmooth<0.001) {
	tf2 = tf*tf2/0.001;
      }
      else {
	tf2 = tf*tf2/fSmooth;
      }
      if (tf2<-0.5) tf2 = -0.5;
      if (tf2>0.5) tf2 = 0.5;
      tf2 = 0.5*sin(PI*tf2) + 0.5; // jumping between notes using horizontally-scaled sine segment
      tf2 = tf2 + ti;
      if ( (tf<0.5 && lowersnap) || (tf>=0.5 && uppersnap) ) {
	outpitch = fAmount*tf2 + ((float)1-fAmount)*outpitch;
      }

      // Add in pitch shift
      outpitch = outpitch + fShift;

      // LFO logic
      tf = fLforate*N/(psAutotalent->noverlap*fs);
      if (tf>1) tf=1;
      psAutotalent->lfophase = psAutotalent->lfophase + tf;
      if (psAutotalent->lfophase>1) psAutotalent->lfophase = psAutotalent->lfophase-1;
      lfoval = psAutotalent->lfophase;
      tf = (fLfosymm + 1)/2;
      if (tf<=0 || tf>=1) {
	if (tf<=0) lfoval = 1-lfoval;
      }
      else {
	if (lfoval<=tf) {
	  lfoval = lfoval/tf;
	}
	else {
	  lfoval = 1 - (lfoval-tf)/(1-tf);
	}
      }
      if (fLfoshape>=0) {
	// linear combination of cos and line
	lfoval = (0.5 - 0.5*cos(lfoval*PI))*fLfoshape + lfoval*(1-fLfoshape);
	lfoval = fLfoamp*(lfoval*2 - 1);
      }
      else {
	// smoosh the sine horizontally until it's squarish
	tf = 1 + fLfoshape;
	if (tf<0.001) {
	  lfoval = (lfoval - 0.5)*2/0.001;
	}
	else {
	  lfoval = (lfoval - 0.5)*2/tf;
	}
	if (lfoval>1) lfoval = 1;
	if (lfoval<-1) lfoval = -1;
	lfoval = fLfoamp*sin(lfoval*PI*0.5);
      }
      // add in quantized LFO
      if (iLfoquant>=1) {
	outpitch = outpitch + (int)(numNotes*lfoval + numNotes + 0.5) - numNotes;
      }


      // Convert back from scale notes to semitones
      outpitch = outpitch + iScwarp; // output scale rotate implemented here
      ti = (int)(outpitch/numNotes + 32) - 32;
      tf = outpitch - ti*numNotes;
      ti2 = (int)tf;
      ti3 = ti2 + 1;
      outpitch = iNote2Pitch[ti3%numNotes] - iNote2Pitch[ti2];
      if (ti3>=numNotes) {
	outpitch = outpitch + 12;
      }
      outpitch = outpitch*(tf - ti2) + iNote2Pitch[ti2];
      outpitch = outpitch + 12*ti;
      outpitch = outpitch - (iNote2Pitch[iScwarp] - iNote2Pitch[0]); //more scale rotation here

      // add in unquantized LFO
      if (iLfoquant<=0) {
	outpitch = outpitch + lfoval*2;
      }


      if (outpitch<-36) outpitch = -48;
      if (outpitch>24) outpitch = 24;

      psAutotalent->outpitch = outpitch;

      //  ---- END Modify pitch in all kinds of ways! ----

      // Compute variables for pitch shifter that depend on pitch
      psAutotalent->inphinc = aref*pow(2,inpitch/12)/fs;
      psAutotalent->outphinc = aref*pow(2,outpitch/12)/fs;
      psAutotalent->phincfact = psAutotalent->outphinc/psAutotalent->inphinc;
    }
    // ************************
    // * END Low-Rate Section *
    // ************************



    // *****************
    // * Pitch Shifter *
    // *****************

    // Pitch shifter (kind of like a pitch-synchronous version of Fairbanks' technique)
    //   Note: pitch estimate is naturally N/2 samples old
    psAutotalent->phasein = psAutotalent->phasein + psAutotalent->inphinc;
    psAutotalent->phaseout = psAutotalent->phaseout + psAutotalent->outphinc;

    //   When input phase resets, take a snippet from N/2 samples in the past
    if (psAutotalent->phasein >= 1) {
      psAutotalent->phasein = psAutotalent->phasein - 1;
      ti2 = psAutotalent->cbiwr - N/2;
      for (ti=-N/2; ti<N/2; ti++) {
	psAutotalent->frag[(ti+N)%N] = psAutotalent->cbf[(ti + ti2 + N)%N];
      }
    }

    //   When output phase resets, put a snippet N/2 samples in the future
    if (psAutotalent->phaseout >= 1) {
      psAutotalent->fragsize = psAutotalent->fragsize*2;
      if (psAutotalent->fragsize > N) {
	psAutotalent->fragsize = N;
      }
      psAutotalent->phaseout = psAutotalent->phaseout - 1;
      ti2 = psAutotalent->cbord + N/2;
      ti3 = (long int)(((float)psAutotalent->fragsize) / psAutotalent->phincfact);
      if (ti3>=N/2) {
	ti3 = N/2 - 1;
      }
      for (ti=-ti3/2; ti<(ti3/2); ti++) {
	tf = psAutotalent->hannwindow[(long int)N/2 + ti*(long int)N/ti3];
	// 3rd degree polynomial interpolator - based on eqns from Hal Chamberlin's book
	indd = psAutotalent->phincfact*ti;
	ind1 = (int)indd;
	ind2 = ind1+1;
	ind3 = ind1+2;
	ind0 = ind1-1;
	val0 = psAutotalent->frag[(ind0+N)%N];
	val1 = psAutotalent->frag[(ind1+N)%N];
	val2 = psAutotalent->frag[(ind2+N)%N];
	val3 = psAutotalent->frag[(ind3+N)%N];
	vald = 0;
	vald = vald - (float)0.166666666667 * val0 * (indd - ind1) * (indd - ind2) * (indd - ind3);
	vald = vald + (float)0.5 * val1 * (indd - ind0) * (indd - ind2) * (indd - ind3);
	vald = vald - (float)0.5 * val2 * (indd - ind0) * (indd - ind1) * (indd - ind3);
	vald = vald + (float)0.166666666667 * val3 * (indd - ind0) * (indd - ind1) * (indd - ind2);
	psAutotalent->cbo[(ti + ti2 + N)%N] = psAutotalent->cbo[(ti + ti2 + N)%N] + vald*tf;
      }
      psAutotalent->fragsize = 0;
    }
    psAutotalent->fragsize++;

    //   Get output signal from buffer
    tf = psAutotalent->cbo[psAutotalent->cbord]; // read buffer

    psAutotalent->cbo[psAutotalent->cbord] = 0; // erase for next cycle
    psAutotalent->cbord++; // increment read pointer
    if (psAutotalent->cbord >= N) {
      psAutotalent->cbord = 0;
    }

    // *********************
    // * END Pitch Shifter *
    // *********************

    ti4 = (psAutotalent->cbiwr + 2)%N;
    if (iFcorr>=1) {
      // The second part of the formant corrector
      // This is a post-filter that re-applies the formants, designed
      //   to result in the exact original signal when no pitch
      //   manipulation is performed.
      // tf is signal input
      // gotta run it 3 times because of a pesky delay free loop
      //  first time: compute 0-response
      tf2 = tf;
      fa = 0;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-psAutotalent->frc[ti])*frlamb + psAutotalent->frb[ti];
	tf = psAutotalent->fbuff[ti][ti4];
	fb = fc - tf*fa;
	psAutotalent->ftvec[ti] = tf*fc;
	fa = fa - psAutotalent->ftvec[ti];
      }
      tf = -fa;
      for (ti=ford-1; ti>=0; ti--) {
	tf = tf + psAutotalent->ftvec[ti];
      }
      f0resp = tf;
      //  second time: compute 1-response
      fa = 1;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-psAutotalent->frc[ti])*frlamb + psAutotalent->frb[ti];
	tf = psAutotalent->fbuff[ti][ti4];
	fb = fc - tf*fa;
	psAutotalent->ftvec[ti] = tf*fc;
	fa = fa - psAutotalent->ftvec[ti];
      }
      tf = -fa;
      for (ti=ford-1; ti>=0; ti--) {
	tf = tf + psAutotalent->ftvec[ti];
      }
      f1resp = tf;
      //  now solve equations for output, based on 0-response and 1-response
      tf = (float)2*tf2;
      tf2 = tf;
      tf = ((float)1 - f1resp + f0resp);
      if (tf!=0) {
	tf2 = (tf2 + f0resp) / tf;
      }
      else {
	tf2 = 0;
      }
      //  third time: update delay registers
      fa = tf2;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-psAutotalent->frc[ti])*frlamb + psAutotalent->frb[ti];
	psAutotalent->frc[ti] = fc;
	psAutotalent->frb[ti] = fb;
	tf = psAutotalent->fbuff[ti][ti4];
	fb = fc - tf*fa;
	fa = fa - tf*fc;
      }
      tf = tf2;
      tf = tf + flpa*psAutotalent->flp;  // lowpass post-emphasis filter
      psAutotalent->flp = tf;
      // Bring up the gain slowly when formant correction goes from disabled
      // to enabled, while things stabilize.
      if (psAutotalent->fmute>0.5) {
	tf = tf*(psAutotalent->fmute - 0.5)*2;
      }
      else {
	tf = 0;
      }
      tf2 = psAutotalent->fmutealph;
      psAutotalent->fmute = (1-tf2) + tf2*psAutotalent->fmute;
      // now tf is signal output
      // ...and we're done messing with formants
    }
    else {
      psAutotalent->fmute = 0;
    }

    // Write audio to output of plugin
    // Mix (blend between original (delayed) =0 and processed =1)
    *(pfOutput++) = (LADSPA_Data) fMix*tf + (1-fMix)*psAutotalent->cbi[ti4];

  }

  // Tell the host the algorithm latency
  *(psAutotalent->m_pfLatency) = (LADSPA_Data) (N-1);
