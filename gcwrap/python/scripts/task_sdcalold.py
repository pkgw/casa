from taskinit import casalog

import asap as sd
from asap._asap import Scantable
from asap.scantable import is_scantable
import sdutil

@sdutil.asaptask_decorator
def sdcalold(infile, antenna, fluxunit, telescopeparam, field, spw, scan, pol, calmode, fraction, noff, width, elongated, markonly, plotpointings, tau, verify, outfile, outform, overwrite, plotlevel):
    with sdutil.sdtask_manager(sdcal_worker, locals()) as worker:
        worker.initialize()
        worker.execute()
        worker.finalize()


class sdcal_worker(sdutil.sdtask_template):
    def __init__(self, **kwargs):
        super(sdcal_worker,self).__init__(**kwargs)
        self.suffix = '_cal'
        # Nothing to be done when calmode='none' and tau=0.0
        if self.calmode=='none' and self.tau==0.0:
            raise Exception("No operation to be done for calmode='none' and tau=0.0. Exiting task.")

    def initialize_scan(self):
        sorg=sd.scantable(self.infile,average=False,antenna=self.antenna)

        if not isinstance(sorg,Scantable):
            raise Exception('Scantable data %s, is not found')

        # A scantable selection
        sel = self.get_selector(sorg)
        sorg.set_selection(sel)
        self.assert_no_channel_selection_in_spw('warn')

        # Copy scantable when usign disk storage not to modify
        # the original table.
        if is_scantable(self.infile) and self.is_disk_storage:
            self.scan = sorg.copy()
        else:
            self.scan = sorg
        del sorg

    def execute(self):
        engine = sdcal_engine(self)
        engine.initialize()

        # apply inputs to scan
        self.set_to_scan()

        # Actual implementation is defined outside the class
        # since those are used in task_sdreduce.
        engine.execute()

        # do opacity (atmospheric optical depth) correction
        sdutil.doopacity(self.scan, self.tau)

        engine.finalize()

    def save(self):
        sdutil.save(self.scan, self.project, self.outform, self.overwrite)


class sdcal_engine(sdutil.sdtask_engine):
    def __init__(self, worker):
        super(sdcal_engine,self).__init__(worker)

    def initialize(self):
        if ( abs(self.plotlevel) > 1 ):
            casalog.post( "Initial Raw Scantable:" )
            self.worker.scan._summary()

    def execute(self):
        scanns = self.worker.scan.getscannos()
        sn=list(scanns)
        casalog.post( "Number of scans to be processed: %d" % (len(sn)) )
        if self.calmode == 'otf' or self.calmode=='otfraster':
            self.__mark()
            if not self.markonly:
                self.worker.scan = sd.asapmath.calibrate( self.worker.scan,
                                                          scannos=sn,
                                                          calmode='ps',
                                                          verify=self.verify )
        else:
            self.worker.scan = sd.asapmath.calibrate( self.worker.scan,
                                                      scannos=sn,
                                                      calmode=self.calmode,
                                                      verify=self.verify )

    def finalize(self):
        if ( abs(self.plotlevel) > 1 ):
            casalog.post( "Final Calibrated Scantable:" )
            self.worker.scan._summary()

        # Plot final spectrum
        if ( abs(self.plotlevel) > 0 ):
            pltfile = self.project + '_calspec.eps'
            sdutil.plot_scantable(self.worker.scan, pltfile, self.plotlevel)

    def __mark(self):
        israster = (self.calmode == 'otfraster')
        marker = sd.edgemarker(israster=israster)
        marker.setdata(self.worker.scan)
        self.npts = self.noff
        marker.setoption(**self.__dict__)
        marker.mark()
        if self.plotpointings:
            marker.plot()
        self.worker.scan = marker.getresult()

