"""
Test sounds and plot waveforms.

This script tests the sound waveform generator for a variety of sounds

"""
import numpy as np
import pyqtgraph as pg
from cnmodel.util import sound
from collections import OrderedDict
import scipy.signal
import sys
import sounddevice as sd
import threading
# try:
#     import PySounds
#     HAVE_PYSOUNDS = True
# except ImportError:
#     HAVE_PYSOUNDS = False


print(sd.query_devices())
class PlaySound():
    def __init__(self):
        pass
        
    def play(self):
        if len(sys.argv) >= 2:
            stimarg = sys.argv[1]
        else:
            exit()

        self.plots = True

        # if HAVE_PYSOUNDS:
        #     PS = PySounds.PySounds()
    
        cf = 2e3
        Fs = 44100  # sample frequency
        level = 80.
        seed = 34978
        fmod = 20.
        dmod = 20.

        if self.plots:
            # waveforms
            win = pg.GraphicsWindow()
            pipwin = win.addPlot(title='sound pip', row=0, col=0)
            pipmodwin = win.addPlot(title='100 \% SAM modulated pip', row=1, col=0)
            noisewin = win.addPlot(title='WB noise', row=2, col=0)
            noisemodwin = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=0)
            clickwin = win.addPlot(title='clicks', row=4, col=0)
            poisson_clickwin = win.addPlot(title='poisson_clicks', row=5, col=0)
            fmwin = win.addPlot(title='fmsweep', row=6, col=0)
            wavewin = win.addPlot(title='wavefile', row=7, col=0)
            # spectra
            pipwins = win.addPlot(title='sound pip Spec', row=0, col=1)
            pipmodwins = win.addPlot(title='100 \% SAM modulated pip', row=1, col=1)
            noisewins = win.addPlot(title='WB noise', row=2, col=1)
            noisemodwins = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=1)
            clickwins = win.addPlot(title='click spec', row=4, col=1)
            poisson_clickwins = win.addPlot(title='poisson_clicks', row=5, col=1)
            fmwins = win.addPlot(title='fmsweep spec', row=6, col=1)
            wavewins = win.addPlot(title='wavefile', row=7, col=1)
        else:
            pipwin = None
            pipmodwin = None
            noisewin = None
            noisemodwin = None
            clickwin = None
            p_clickwin = None
            pipwins = None
            pipmodwins = None
            noisewins = None
            noisemodwins = None
            clickwins = None
            p_clickwins = None
            fmwins = None

        stims = OrderedDict(
            [
                ("pip", (pipwin, sound.TonePip)),
                ("pipmod", (pipmodwin, sound.SAMTone)),
                ("noise", (noisewin, sound.NoisePip)),
                ("noisemod", (noisemodwin, sound.SAMNoise)),
                ("clicks", (clickwin, sound.ClickTrain)),
                ("poisson_clicks", (poisson_clickwin, sound.ClickTrain)),
                ("wavefile", (wavewin, sound.ReadWavefile)),
            ]
        )

        specs = OrderedDict(
            [
                ("pip", (pipwins, sound.TonePip)),
                ("pipmod", (pipmodwins, sound.SAMTone)),
                ("noise", (noisewins, sound.NoisePip)),
                ("noisemod", (noisemodwins, sound.SAMNoise)),
                ("clicks", (clickwins, sound.ClickTrain)),
                ("poisson_clicks", (poisson_clickwins, sound.ClickTrain)),
                ("wavefile", (wavewins, sound.ReadWavefile)),
            ]
        )

        if stimarg == 'all':
            stimlist = list(stims.keys())
        else:
            stimlist = [stimarg]
        for stim in stimlist:
            print(stim)
            if stim in ['clicks']:
                self.wave = stims[stim][1](rate=Fs, duration=1.0, dbspl=level,
                                 click_duration=1e-4, click_starts=1e-3*np.linspace(10, 500, 10))
            elif stim in ["poisson_clicks"]:
                clickrate = 50.0 # Hz
                eventintervals = np.random.exponential(
                        1.0 / clickrate, int(1.0*clickrate)
                    )

                events = np.cumsum(eventintervals)
                events = events[events < 0.8]
                self.wave = stims[stim][1](
                    rate=Fs,
                    duration=1.0,
                    dbspl=level,
                    click_duration=1e-4,
                    click_starts=events,
                )
                # wave = stims[
            elif stim in ['fmsweep']:
                self.wave = stims[stim][1](rate=Fs, duration=0.5, dbspl=level,
                                    start=0., ramp='linear', freqs=[16000, 200])
            elif stim in ['pip', 'pipmod', 'noise', 'noisemod']:
                self.wave = stims[stim][1](rate=Fs, duration=2.0, f0=cf, dbspl=level, 
                                 pip_duration=1.8, pip_start=[10e-3], ramp_duration=2.5e-3,
                                 fmod=fmod, dmod=dmod, seed=seed)
            if self.plots:
                stims[stim][0].plot(self.wave.time, self.wave.sound)
                f, Pxx_spec = scipy.signal.periodogram(self.wave.sound, Fs) #, window='flattop', nperseg=8192,
                                   # noverlap=512, scaling='spectrum')
                specs[stim][0].plot(f, np.sqrt(Pxx_spec))
        
            event = threading.Event()
            samplerate = sd.query_devices(1, 'output')['default_samplerate']
            print("samplerate: ", samplerate)
            self.current_frame = 0
            # def callback(self, outdata, frames, time, status):
            #         global current_frame
            #         if status:
            #             print(status)
            #         data = wave.sound
            #         chunksize = min(len(data) - current_frame, frames)
            #         outdata[:chunksize] = data[current_frame:current_frame + chunksize]
            #         if chunksize < frames:
            #             outdata[chunksize:] = 0
            #             raise sd.CallbackStop()
            #         current_frame += chunksize
            stream = sd.OutputStream(device=1, channels=1, callback=self.callback,
                                 samplerate=samplerate)
            with stream:
                event.wait(timeout=1.5)
            print("done playing")
        # if self.plots and sys.flags.interactive == 0:
        #      pg.QtWidgets.QApplication.exec()
            # if HAVE_PYSOUNDS:
            #
            #     print ('Playing %s' % stim)
            #
            #     PS.playSound(wave.sound, wave.sound, Fs)
    def callback(self, outdata, frames, time, status):
            # global current_frame
            if status:
                print(status)
            data = np.array(self.wave.sound)
            chunksize = min(len(data) - self.current_frame, frames)
            outdata[:chunksize] = np.array(data[self.current_frame:self.current_frame + chunksize])[...,np.newaxis]
            if chunksize < frames:
                outdata[chunksize:] = 0
                raise sd.CallbackStop()
            self.current_frame += chunksize

         
if __name__ == '__main__':
    # if not HAVE_PYSOUNDS:
    #     print("Could not import PySounds; will not play audio.")
    P = PlaySound()
    P.play()

    # if P.plots : # and sys.flags.interactive == 0:
    #      pg.QtWidgets.QApplication.exec()