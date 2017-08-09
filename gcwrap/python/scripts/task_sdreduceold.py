from taskinit import casalog

import asap as sd
import sdutil

import task_sdcalold as task_sdcalold
import task_sdaverageold as task_sdaverageold
import task_sdbaselineold as task_sdbaselineold

@sdutil.asaptask_decorator
def sdreduceold(infile, antenna, fluxunit, telescopeparam, field, spw, restfreq, frame, doppler, timerange, scan, pol, calmode, fraction, noff, width, elongated, markonly, plotpointings, tau, average, timeaverage, tweight, scanaverage, averageall, polaverage, pweight, kernel, kwidth, chanwidth, maskmode, thresh, avg_limit, edge, blfunc, order, npiece, applyfft, fftmethod, fftthresh, addwn, rejwn, clipthresh, clipniter, verifycal, verifysm, verifybl, verbosebl, bloutput, blformat, showprogress, minnrow, outfile, outform, overwrite, plotlevel):
    with sdutil.sdtask_manager(sdreduce_worker, locals()) as worker:
        worker.initialize()
        worker.execute()
        worker.finalize()


class sdreduce_worker(sdutil.sdtask_template):
    def __init__(self, **kwargs):
        super(sdreduce_worker,self).__init__(**kwargs)
        self.suffix = '_cal'

    def initialize_scan(self):
        # instantiate scantable
        self.scan = sd.scantable(self.infile, average=False, antenna=self.antenna)

        # restorer
        self.restorer = sdutil.scantable_restore_factory(self.scan,
                                                         self.infile,
                                                         self.fluxunit,
                                                         '', # specunit=''
                                                         self.frame,
                                                         self.doppler,
                                                         self.restfreq)

        # Apply selection
        self.scan.set_selection(self.get_selector())

    def execute(self):
        # calibration stage
        casalog.post( "*** sdcalold stage ***" )
        self.verify = self.verifycal
        engine = task_sdcalold.sdcal_engine(self)
        engine.initialize()
        engine.execute()
##         self.scan = engine.get_result()
##         task_sdcalold.prior_plot(self.scan, self.plotlevel)
##         self.scan = task_sdcalold.docalibration(self.scan, self.calmode,
##                                              self.fraction,
##                                              self.noff, self.width,
##                                              self.elongated,
##                                              self.markonly, self.plotpointings,
##                                              self.verifycal)

        # apply input parameters
        self.set_to_scan()

        # opacity correction
        sdutil.doopacity(self.scan, self.tau)

        ## channel splitting
        #sdutil.dochannelrange(self.scan, self.channelrange)

        #WORKAROUND for new tasks (in future this should be done in sdutil)
        if not self.timeaverage: self.scanaverage = False
        # averaging stage
        if self.average:
            self.scan = sdutil.doaverage(self.scan, self.scanaverage,
                                         self.timeaverage, self.tweight,
                                         self.polaverage, self.pweight,
                                         self.averageall)
        else:
            casalog.post( "No averaging was applied..." )
##         task_sdcalold.posterior_plot(self.scan, self.project, self.plotlevel)
        engine.finalize()
        del engine

        # smoothing stage
        casalog.post( "" )
        casalog.post( "*** sdsmooth stage ***" )
        if self.kernel.lower() not in ['none', '']:
            self.verify = self.verifysm
            engine = task_sdaverageold.sdsmooth_engine(self)
            engine.initialize()
            engine.execute()
##             self.scan = engine.get_result()
            engine.finalize()
            del engine
        else:
            casalog.post( "No smoothing was applied..." )


        # baseline stage
        casalog.post( "" )
        casalog.post( "*** sdbaselineold stage ***")
        if self.blfunc != 'none':
            self.verify = self.verifybl
            self.verbose = self.verbosebl
            engine = task_sdbaselineold.sdbaseline_engine(self)
            engine.initialize()
            engine.execute()
##             self.scan = engine.get_result()
            engine.finalize()
            del engine
        else:
            casalog.post( "No baseline subtraction was applied..." )

    def save(self):
        # write result on disk
        sdutil.save(self.scan, self.project, self.outform, self.overwrite)

    def cleanup(self):
        # restore
        if hasattr(self,'restorer') and self.restorer:
            self.restorer.restore()
