
this package contains the software used by the Sonification of Solar 
Harmonics (SoSH) Project.  to run it, you will need to install pure 
data, which can be found at http://puredata.info/downloads/pure-data or 
at http://msp.ucsd.edu/software.html .  once installed, make sure it is 
working with your soundcard by clicking on "Test Audio and MIDI" under 
the "Media" menu.

NOTE: pure data does not work in windows with block sizes above 65536, 
whereas the block size used here may be as large 262144.  this has 
worked fine for mac and linux.  when running in windows, make sure the 
block size is set appropriately.  a discussion and possible fix of this 
issue can be found at https://github.com/pure-data/pure-data/issues/68 .

files to note:
modefilter_standalone.pd - the stand alone version of the patch.
modefilter.pd - version to be used as a subpatch.
instructions_audio.txt - an extensive introduction and explanation of 
the patches.
quickstart_audio.txt - much shorter instructions.
example_sum.pd - patch to play 5 tones at once or in sequence.
example_sequencer.pd - patch to combine arbitrary lists of modes and 
concatenate them.

files written by the patches will go into the "wav_out" subdirectory.

all included files:
applywindow.pd        fft-resynth-negm.pd       qlist.hmi
arbitrarySRmulti.pd   fft-resynth-posm.pd       qlist.mdi
arbitrarySR.pd        instructions_audio.txt    qlist.test8
audio_safety~.pd      loadaudio.pd              quickstart_audio.txt
calcbinshift.pd       makegain.pd               README_audio.txt
example1.pd           makeoutputstring.pd       set-directory.pd
example2.pd           modeaddition.pd           symbolchange.pd
example3.pd           modecat.pd                text-file-reader.pd
example_addition.pd   modefilter0.pd            triggerlogic.pd
example_concat2.pd    modefilter.pd             vuzi.pd
example_concat.pd     modefilter_standalone.pd  window-gen.pd
example_sequencer.pd  modesum.pd                
example_sum.pd        numberchange.pd
fft-analysis.pd       parsedaylnm.pd
