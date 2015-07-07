import datetime
import sys
import time

class ProgressBar:
    def __init__(self, report_frequency):
        self.start_time = time.time()
        self.last_report_time = self.start_time
        self.report_frequency = report_frequency

    def GetEta(self, progress):
        elapsed = time.time() - self.start_time
        if progress > 0:
            return elapsed / progress - elapsed
        else:
            return 0

    def GetEtaAsString(self, progress):
        eta = self.GetEta(progress)
        t = datetime.datetime.now() + datetime.timedelta(seconds=eta)
        return t.strftime('%Y-%m-%d %H:%M:%S')

    def MaybeReport(self, progress):
        if self.report_frequency == 0:
            return
        if time.time() - self.last_report_time > self.report_frequency:
            percent = '%.2f %%' % (progress * 100)
            sys.stderr.write(percent + ' ' +
                             self.GetEtaAsString(progress) + '\n')
            self.last_report_time = time.time()
