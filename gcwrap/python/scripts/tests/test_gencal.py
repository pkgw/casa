import os
import sys
import shutil
import subprocess
import numpy
import numpy.ma as ma
import testhelper as th
from __main__ import default
from tasks import gencal
from taskinit import *
import unittest

'''
Unit tests for gencal 
'''
#
# ToDo:
# add more tests
# once more independent tests (e.g. comparison
# the AIPS REWAY results) add reference mses
# and do tests against them
# 

datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/gencal/'

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/gencal/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('gencal tests will use data from '+datapath)         


class gencal_antpostest(unittest.TestCase):

    # Input and output names
    msfile = 'tdem0003gencal.ms'
#    if testmms:
#        msfile = 'tdem0003gencal.mms'
    caltable = 'anpos.cal'
    reffile1 = datapath+'anpos.manual.cal'
    reffile2 = datapath+'anpos.auto.cal'
    res = False

    def setUp(self):
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)

        shutil.copytree(datapath+self.msfile, self.msfile, symlinks=True)

    def tearDown(self):
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)

        shutil.rmtree(self.caltable,ignore_errors=True)

    def test_antpos_manual(self):
        """
        gencal: test manual antenna position correction 
        """
        gencal(vis=self.msfile,
               caltable=self.caltable, 
               caltype='antpos',
               antenna='ea12,ea22',
               parameter=[-0.0072,0.0045,-0.0017, -0.0220,0.0040,-0.0190])

        self.assertTrue(os.path.exists(self.caltable))

        # ToDo:check generated caltable. Wait for new caltable
        
        # Compare with reference file from the repository
        reference = self.reffile1
        self.assertTrue(th.compTables(self.caltable, reference, ['WEIGHT','OBSERVATION_ID']))

    def test_antpos_auto_evla(self):
        """
        gencal: test automated antenna position correction
        """
        # check if the URL is reachable
        import urllib.request, urllib.error, urllib.parse
        # current EVLA baseline correction URL
        evlabslncorrURL="http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year="
        try: 
          urlaccess=urllib.request.urlopen(evlabslncorrURL+"2010", timeout=30.0) 
          gencal(vis=self.msfile,
                 caltable=self.caltable,
                 caltype='antpos',
                 antenna='',
                 parameter='')

          self.assertTrue(os.path.exists(self.caltable))
          
          # ToDo: check for generated caltable
          
          # Compare with reference file from the repository
          reference = self.reffile2
          self.assertTrue(th.compTables(self.caltable, reference, ['WEIGHT','OBSERVATION_ID']))

        except urllib.error.URLError as err:
          print("Cannot access %s , skip this test" % evlabslncorrURL)
          self.res=True


class test_gencal_antpos_alma(unittest.TestCase):
    """
    Tests the automatic generation of antenna position corrections for ALMA

    New REST web service:
    https://bitbucket.sco.alma.cl/projects/ALMA/repos/almasw/browse/CONTROL-SERVICES/PositionsService


    Old SOAP web service:
    http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl
    Example minimalistic use of a client to query the service:
      from suds.client import Client
      srv_wsdl_url = 'http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl'
      ws_cli = Client(srv_wsdl_url)
      resp = ws_cli.service.getAntennaPositions("CURRENT.AOS", "DA49",
                                                "2017-01-30T01:53:54")
    """

    # setup of the ALMA TMC DB AntennaPadService
    ALMA_SRV_WSDL_URL = 'http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl'

    # For this MS, there is position information for 25 out of the 29 antennas
    # (at 2013-11-15T10:26:19)
    ALMA_MS = 'uid___A002_X72c4aa_X8f5_scan21_spw18_field2_corrXX.ms'
    CAL_TYPE = 'antpos'
    REF_CALTABLE_MANUAL = os.path.join(datapath,
                                       'alma_reference/A002_X72c4aa_ref_ant_pos.manual.cal')
    REF_CALTABLE_AUTO = os.path.join(datapath,
                                     'alma_reference/A002_X72c4aa_ref_ant_pos.auto.cal')
    IGNORE_COLS = ['WEIGHT','OBSERVATION_ID']

    def setUp(self):
        if (os.path.exists(self.ALMA_MS)):
            shutil.rmtree(self.ALMA_MS)

        flagdata_datapath = os.path.join(datapath, '../flagdata/')
        shutil.copytree(os.path.join(flagdata_datapath, self.ALMA_MS),
                        self.ALMA_MS, symlinks=True)

    def tearDown(self):
        if (os.path.exists(self.ALMA_MS)):
            shutil.rmtree(self.ALMA_MS)

    def remove_caltable(self, ct_name):
        """ Removes a cal table. ct_name: path to the caltable """
        import shutil
        shutil.rmtree(ct_name)

    def test_antpos_alma_manual(self):
        """
        gencal: manual antenna position correction on ALMA table
        """

        out_caltable = 'ant_pos_man.cal'
        gencal(vis=self.ALMA_MS,
               caltable=out_caltable,
               caltype=self.CAL_TYPE,
               antenna='DV07,DV10,DV11',
               parameter=[-0.0072,0.0045,-0.0017, -0.0220,0.0040,-0.0190])

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created")

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_MANUAL,
                                      self.IGNORE_COLS))

        self.remove_caltable(out_caltable)

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_alma_server_SOAP_methods(self):
        """
        gencal: connection to alma TCM DB AntennaPadService for ALMA
        """
        try:
            import urllib.request, urllib.error, urllib.parse
            from suds.client import Client
            ws_cli = Client(self.ALMA_SRV_WSDL_URL)

            # Basic check that the schema has the minimum requirement
            method_name = 'getAntennaPositions'
            self.assertTrue(callable(getattr(ws_cli.service, method_name)),
                            'The client service should have this method: {}, and '
                            'it should be callable.'.format(method_name))
        except ImportError as exc:
            print('Cannot import required dependencies to query the ALMA TCM DB '
                  'web service')
            raise
        except urllib.error.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web'
                  'service')
            raise

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_auto_alma_SOAP_empty_query(self):
        """
        gencal: empty query (empty antennas list) to the (old) SOAP TCMDB AntennaPadService
        web service (ALMA)
        """
        try:
            import correct_ant_posns_alma as almacor

            resp = almacor.query_tmcdb_antennas_rest([], '2017-01-01T16:53:54.000')
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB '
                  'web service')
            raise
        except urllib.error.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web'
                  'service')
            raise

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_auto_web_srv_SOAP_alma(self):
        """
        gencal: auto gencal using data from TCM DB AntennaPadService (ALMA)
        """

        import urllib.request, urllib.error, urllib.parse

        out_caltable = 'ant_pos_web_srv.cal'
        try:
            # This will import the required libraries, urllib2, suds, etc.
            # Coul also use additional parameters: antenna='', parameter=''
            gencal(vis=self.ALMA_MS, caltable=out_caltable, caltype=self.CAL_TYPE)
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB '
                  'web service')
            raise
        except urllib.error.URLError:
            print('Connection/network error while querying the ALMA TCM DB web'
                  'service')
            raise

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created: {0}".
                        format(out_caltable))

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_AUTO,
                                      self.IGNORE_COLS))
        self.remove_caltable(out_caltable)


    @unittest.skip('REST Position service needs validation and final deployment')
    def tmp_disabled_test_antpos_auto_alma_REST_empty_query(self):
        """
        gencal: empty query (empty antennas list) to the (new) REST TCMDB Positions
        web service (ALMA)
        """
        import urllib.request, urllib.error, urllib.parse

        TEST_HOSTNAME = 'https://2018may.asa-test.alma.cl'

        hostname = TEST_HOSTNAME
        port = 443
        api = 'antenna-position/position/antenna'
        try:
            import requests
            import correct_ant_posns_alma as almacor

            tstamp = '2017-01-01T16:53:54.000'
            # query via correct_ant_posns function
            resp = almacor.query_tmcdb_antennas_rest([], tstamp)
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))

            # query directly via requests
            url = '{}:{}/{}?antenna={}&timestamp={}'.format(hostname, port, api, '',
                                                            '2017-01-01T16:53:54.000')
            resp = requests.get(url)
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB '
                  'web service')
            raise
        except urllib.error.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web'
                  'service')
            raise

    @unittest.skip('REST Position service needs validation and final deployment')
    def tmp_disabled_test_antpos_auto_web_srv_REST_alma(self):
        """
        gencal: auto gencal using data from TCMDB Positions service (ALMA)
        """

        import urllib.request, urllib.error, urllib.parse

        out_caltable = 'ant_pos_web_srv.cal'
        try:
            # This will import the required libraries, urllib2, suds, etc.
            # Coul also use additional parameters: antenna='', parameter=''
            gencal(vis=self.ALMA_MS, caltable=out_caltable, caltype=self.CAL_TYPE)
        except urllib.error.URLError:
            print('Connection/network error while querying the ALMA TCMDB Positions web'
                  'service')
            raise

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created: {0}".
                        format(out_caltable))

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_AUTO,
                                      self.IGNORE_COLS))
        self.remove_caltable(out_caltable)


def suite():
    return [gencal_antpostest,
            test_gencal_antpos_alma]
