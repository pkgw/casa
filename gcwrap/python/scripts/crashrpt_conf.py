import os
import tempfile
import getpass

temporaryDirectory = tempfile.gettempdir() + "/casa-crashreporter/" + getpass.getuser() + "/" + str(os.getpid())
