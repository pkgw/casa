import os
import time
import threading
from casac import casac
class TelemetryLogMonitor:

    def __init__(self):
        self.showWarning = True

    # Limit in kilobytes
    def isWithinLimit(self, filename, limit, casa):
        try:
            filesize = os.path.getsize(filename)/1024
            if (filesize > limit):
                if (self.showWarning):
                    print("Logfile size is too large. Disabling telemetry.")
                    print("Filesize: " + str(filesize))
                    print("Limit: " + str(limit))
                casa['state']['telemetry-enabled'] = False
                self.showWarning = False
        except OSError:
            # If file doesn't exist, we'll create a new one
            pass
        # Else if telemetry disabled enable it again?

    def monitorLogFileSize(self, filename, limit, interval, casa):
        while (True):
            self.isWithinLimit (filename, limit, casa)
            time.sleep(interval)

    def start(self, filename, limit, interval, casa):
        try:
            logmonitor_thread = threading.Thread(target=self.monitorLogFileSize, args=[filename, limit, interval, casa])
            logmonitor_thread.start()
        except (KeyboardInterrupt, SystemExit):
            cleanup_stop_thread()
            sys.exit()
