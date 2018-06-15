import os
import fnmatch
import subprocess
import tarfile
from casac import casac
import datetime
import time
import urllib2

class telemetry:

    def __init__(self, logdir, logpattern, casa):
        if (logdir is None or logdir == ""):
            print "No log directory provided. Can't submit telemetry data."
            return
        if (logpattern is None or logpattern == ""):
            print "No log pattern provided. Can't submit telemetry data."
            return
        self.logdir = logdir
        self.logpattern = logpattern
        self.stampfile = self.logdir + "/telemetry.stamp"
        self.casa = casa


    def setCasaLog(self, logger):
        self.logger = logger

    def submitStatistics(self):
        if (self.casa['state']['telemetry-enabled'] == True):
            self.logger.poststat("Stop CASA")
            self.logger.post("Checking telemetry submission interval")
            self.createStampFile()
            if (self.isSubmitInterval()):
                self.send('https://casa.nrao.edu/cgi-bin/crash-report.pl')
                self.refreshStampFile()


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
        try:
            urllib2.urlopen('https://casa.nrao.edu/', timeout=2)
        except urllib2.URLError as err:
            return

        # Find logfiles
        for file in os.listdir(self.logdir):
            if fnmatch.fnmatch(file, self.logpattern):
                #print "Matched: " + file
                logfiles.append(file)

        if (len(logfiles) > 0):
            #Tar logfiles
            current_date = datetime.datetime.today().strftime('%Y%m%d%H%M%S')
            tarfileid = self.logdir + "telemetry-" \
                        + telemetryhelper.getUniqueId() + "-" \
                        + current_date + ".tar.gz"
            tar = tarfile.open(tarfileid, "w:gz")

            for logfile in logfiles:
                tar.add(self.logdir + "/" + logfile,
                        arcname='telemetry/'+logfile)
            tar.close()

            file_param = 'file=@' + tarfileid #+ '\"'
            # Submit tarfile
            #print ['curl', '-F', file_param , telemetry_url]
            proc = subprocess.Popen(['curl', '-F', file_param , telemetry_url],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            cmd_out, cmd_err = proc.communicate()
            if cmd_out != None:
                self.logger.post(cmd_out, 'DEBUG1')
            if cmd_err != None:
                self.logger.post(cmd_err, 'DEBUG1')

            # Remove files
            for logfile in logfiles:
                os.remove(self.logdir + "/" + logfile)
                #print "Removed " + self.logdir + "/" + logfile
            os.remove(tarfileid)
            self.logger.post("Removed" + tarfileid)
        else:
            self.logger.post("No telemetry files to submit.")
