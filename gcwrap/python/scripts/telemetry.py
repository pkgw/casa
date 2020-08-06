import os
import fnmatch
import subprocess
import tarfile
from casac import casac
import datetime
import time
import urllib.request, urllib.error, urllib.parse
import __casac__
import TelemetryLogMonitor
import ssl

class telemetry:

    def __init__(self, casa):

        self.setCasaVersion()
        self.setHostId()
        self.logdir = casa['dirs']['rc']
        self.variantSuffix = ""
        if len(casa['variant'])>1:
            self.variantSuffix = "-" + casa['variant']
        self.logpattern = 'casastats-' + self.casaver + '-' + self.hostid + '*' + self.variantSuffix + '.log'
        self.sendlogpattern = 'casastats-*'+ self.hostid + '*.log'
        self.stampfile = self.logdir + '/telemetry-' + self.hostid + '.stamp'
        self.casa = casa

        logfiles = []

        # Check if user has defined a telemetry log location
        casa_util = __casac__.utils.utils()

        if (casa_util.getrc("TelemetryLogDirectory") != 'Unknown value'):
            self.logdir = casa_util.getrc("TelemetryLogDirectory")

        for file in os.listdir(self.logdir):
            if fnmatch.fnmatch(file, self.logpattern):
                 #print "Matched: " + file
                 logfiles.append(file)

        logfiles.sort(reverse=True)
        # Size of the existing (non-active) logfiles
        inactiveTLogSize = 0

        if (logfiles and logfiles[0] != None):
            print("Found an existing telemetry logfile: " + self.logdir  + "/" + logfiles[0])
            casa['files']['telemetry-logfile'] = self.logdir  + "/" + logfiles[0]
            for i in range(1, len(logfiles)):
                inactiveTLogSize = inactiveTLogSize + os.path.getsize(self.logdir  + "/" + logfiles[i])/1024
                #print "Inactive log size: " + str(inactiveTLogSize)
        else :
            print("Creating a new telemetry file")
            self.setNewTelemetryFile()

        # Setup Telemetry log size monitoring
        # Size limit for the telemetry logs
        tLogSizeLimit = 20000
        # File size check interval
        tLogSizeInterval = 60
        try:
            tLogSizeLimit = int(casa_util.getrc("TelemetryLogLimit"))
            tLogSizeInterval = int(casa_util.getrc("TelemetryLogSizeInterval"))
        except:
            pass
        # Subtract the inactive log sizes from the total log file size limit
        tLogSizeLimit = tLogSizeLimit - inactiveTLogSize
        if (tLogSizeLimit <= 0):
            print("Telemetry log size limit exceeded. Disabling telemetry.")
            casa['state']['telemetry-enabled'] = False
        else :
            tLogMonitor = TelemetryLogMonitor.TelemetryLogMonitor()
            tLogMonitor.start(casa['files']['telemetry-logfile'],tLogSizeLimit, tLogSizeInterval, casa)
            print("Telemetry initialized. Telemetry will send anonymized usage statistics to NRAO.")
            print('You can disable telemetry by adding the following line to your ~/.casarc file:')
            print('EnableTelemetry: False')

    def setNewTelemetryFile(self):
        self.casa['files']['telemetry-logfile'] =  self.logdir + '/casastats-' + self.casaver +'-'  + self.hostid + "-" + time.strftime("%Y%m%d-%H%M%S", time.gmtime()) + self.variantSuffix + '.log'
        # Work around the chicken/egg problem with telemetry/logger initialization
        if hasattr(self, 'logger'):
            self.logger.setstatslogfile(self.casa['files']['telemetry-logfile'])

    def setCasaVersion(self):
        myUtils = casac.utils()
        ver = myUtils.version()
        self.casaver = str(ver[0])+ str(ver[1]) + str(ver[2])+ "-" + str(ver[3])

    def setHostId(self):
        telemetryhelper = casac.telemetryhelper()
        self.hostid = telemetryhelper.getUniqueId()

    def setCasaLog(self, logger):
        self.logger = logger

    def submitStatistics(self):
        if (self.casa['state']['telemetry-enabled'] == True):
            self.logger.post("Checking telemetry submission interval")
            self.createStampFile()
            if (self.isSubmitInterval()):
                postingUrl = 'https://casa.nrao.edu/cgi-bin/crash-report.pl'
                if 'CASA_CRASHREPORT_URL' in os.environ :
                    postingUrl = os.environ['CASA_CRASHREPORT_URL']
                self.send(postingUrl)
                self.refreshStampFile()
                self.setNewTelemetryFile()


    def isSubmitInterval(self):
        currentTime = time.time()
        lastUpdateTime = time.time()
        if (os.path.isfile(self.stampfile)):
            lastUpdateTime = os.path.getmtime(self.stampfile)

        # Check update checkSubmitInterval
        interval = 604800
        utils = casac.utils()
        if (utils.getrc("TelemetrySubmitInterval") != 'Unknown value'):
            interval = float(utils.getrc("TelemetrySubmitInterval"))
        if ((currentTime - lastUpdateTime)> interval):
            self.logger.post("Telemetry submit interval reached, submitting telemetry data.")
            return True
        else:
            self.logger.post("Telemetry submit interval not reached. Not submitting data.")
            #print "lastUpdateTime" +str(lastUpdateTime)
            #print "currentTime" +str(currentTime)
            self.logger.post("Next telemetry data submission in: " + str(datetime.timedelta(  \
                    seconds=(interval-(currentTime-lastUpdateTime)))))
            return False

    def createStampFile(self):
        #print "Checking for stampfile " + self.stampfile
        if not os.path.isfile(self.stampfile):
            self.logger.post("Creating a new telemetry time stamp file." + self.stampfile)
            open(self.stampfile, 'a').close()

    def refreshStampFile(self):
        os.utime(self.stampfile, None)

    def send(self, telemetry_url):

        telemetryhelper = casac.telemetryhelper()
        logfiles = []

        # Test if internet connection is available.
        context = ssl._create_unverified_context()
        try:
            urllib.request.urlopen('https://casa.nrao.edu/', timeout=20, context=context)
        except urllib.error.URLError as err:
            self.logger.post("No telemetry server available. Not submitting data")
            return

        # Find logfiles
        for file in os.listdir(self.logdir):
            if fnmatch.fnmatch(file, self.sendlogpattern):
                #print "Matched: " + file
                logfiles.append(file)

        if (len(logfiles) > 0):
            #Tar logfiles
            current_date = datetime.datetime.today().strftime('%Y%m%d%H%M%S')
            tarfileid = self.logdir + "/telemetry-" \
                            + telemetryhelper.getUniqueId() + "-" \
                            + current_date + ".tar.gz"
            try:
                tar = tarfile.open(tarfileid, "w:gz")
                for logfile in logfiles:
                    tar.add(self.logdir + "/" + logfile,
                            arcname='telemetry/'+logfile)
                tar.close()
            except Exception as e:
                self.logger.post("Couldn't create telemetry tarfile")
                self.logger.post(str(e))

            try:
                file_param = 'file=@' + tarfileid #+ '\"'
                # Submit tarfile
                #print ['curl', '-F', file_param , telemetry_url]
                proc = subprocess.Popen(['curl', '-F', file_param , telemetry_url],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                cmd_out, cmd_err = proc.communicate()
                if cmd_out != None:
                    self.logger.post(cmd_out, 'DEBUG1')
                if cmd_err != None:
                    self.logger.post(cmd_err, 'DEBUG1')
            except Exception as e:
                self.logger.post("Couldn't submit telemetry logs")
                self.logger.post(str(e))

            # Remove files
            for logfile in logfiles:
                try:
                    os.remove(self.logdir + "/" + logfile)
                except Exception as e:
                    self.logger.post("Couldn't remove logfile " + self.logdir + "/" + logfile)
                    self.logger.post(str(e))
                    #print "Removed " + self.logdir + "/" + logfile
            try:
                os.remove(tarfileid)
                self.logger.post("Removed" + tarfileid)
            except Exception as e:
                    self.logger.post("Couldn't remove  " + tarfileid)
                    self.logger.post(str(e))
        else:
            self.logger.post("No telemetry files to submit.")
