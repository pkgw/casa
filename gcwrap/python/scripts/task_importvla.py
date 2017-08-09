import os
from taskinit import *
from casa_stack_manip import stack_frame_find
def importvla(archivefiles,vis,
              bandname,frequencytol,
              project,
              starttime,stoptime,
              applytsys,
              autocorr,antnamescheme,keepblanks,evlabands):
        i=0
        overwrite=True
        ok = True
        try:
                casalog.origin('importvla')
                if ((type(vis)!=str) | (vis=='') | (os.path.exists(vis))):
                        raise Exception('Need valid visibility file name (bad name or already exists)')
                if (os.path.exists(vis)): raise Exception('Visibility file already exists - remove or rename')
                for archivefile in archivefiles:
                        if i>0: overwrite=False
                        myf = stack_frame_find( )
                        vftask = myf['vftask']
                        vlafiller = vftask.fill
                        if ((type(archivefile)==str) & (os.path.exists(archivefile))):
                                vlafiller(msname=vis,inputfile=archivefile,
                                          overwrite=overwrite,
                                          bandname=bandname,freqtol=frequencytol,
                                          project=project, start=starttime,
                                          stop=stoptime, applytsys=applytsys,
                                          keepautocorr=autocorr,
                                          antnamescheme=antnamescheme,
                                          keepblanks=keepblanks,
                                          evlabands=evlabands)
                                i += 1
                        else:
                                raise Exception('Archive file not found - please verify the name')
        except Exception as instance:
                print('*** Error importing %s to %s:' % (archivefiles, vis))
                raise Exception(instance)

    # Write history
        try:
                param_names = importvla.__code__.co_varnames[:importvla.__code__.co_argcount]
                param_vals = [eval(p) for p in param_names]
                ok &= write_history(mstool(), vis, 'importvla', param_names,
                                    param_vals, casalog)
        except Exception as instance:
                casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                             'WARN')

    # write initial flag version
        try:
                aflocal = casac.agentflagger()
                ok &= aflocal.open(vis);
                ok &= aflocal.saveflagversion('Original',
                                                        comment='Original flags at import into CASA',
                                                        merge='replace')
                ok &= aflocal.done();
        except Exception as instance:
                print('*** Error writing initial flag version of %s:' % vis)
                raise Exception(instance)
