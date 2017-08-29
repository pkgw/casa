import os
import tempfile
import getpass

systemTempDir = tempfile.gettempdir()
temporaryDirectoryCommon = tempfile.gettempdir() + "/casa-crashreporter/"
temporaryDirectory = temporaryDirectoryCommon + getpass.getuser() + "/" + str(os.getpid())
