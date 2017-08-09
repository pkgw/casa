

class callibrary(object):

    def __init__(self):
        self.cld={}   # as dict

    def __repr__(self):
        return "<callibrary>"

    def cldefinstance(self):
        definst={
            "field" :"",
            "intent":"",
            "spw": "",
            "obs": "",
            "fldmap" : [],
            "obsmap" : [],
            "spwmap" : [],
            "antmap" : [],
            #"calwt" : False,
            "tinterp" : "",
            "finterp" : "",
            "reach" : ""
            }
        return definst


    def clear(self):
        self.__init__()

    def len(self):
        return len(self.cld)

    def addold(self,field='',spw='',intent='',
               gaintable='',gainfield='',interp='',spwmap=[],calwt=False):

        if len(gaintable)<1:
            raise Exception('Please specify at least a gaintable.')

        # insist all cal params are lists
        #  NB: data selection params are _not_ lists
        if (not isinstance(gaintable,list)):
            gaintable=[gaintable]

        if (not isinstance(gainfield,list)):
            gainfield=[gainfield]
        if (not isinstance(interp,list)):
            interp=[interp]
        if (not isinstance(calwt,list)):
            calwt=[calwt]
        if isinstance(spwmap,list) and len(spwmap)>0:
            if (not isinstance(spwmap[0],list)):
                spwmap=[spwmap]  # nest it
        else:
            spwmap=[]


        for itab in range(len(gaintable)):
            tint='linear'
            fint=''
            sinterp=interp[itab].split(',') if itab<len(interp) else []
            if len(sinterp)>0 and len(sinterp[0])>0:
                tint=sinterp[0]
                fint=sinterp[1] if (len(sinterp)>1) else ''
            self.add(caltable=gaintable[itab],
                     field=field,
                     spw=spw,
                     intent=intent,
                     tinterp=tint,
                     finterp=fint,
                     calwt=calwt[itab] if itab<len(calwt) else calwt[len(calwt)-1],
                     fldmap=gainfield[itab] if itab<len(gainfield) else '',
                     spwmap=spwmap[itab] if itab<len(spwmap) else []
                     )
            

    # 
    def add(self,caltable,
            field='',intent='',spw='',obs='',
            tinterp='',finterp='',reach='',calwt=True,
            obsmap=[],fldmap=[],spwmap=[],antmap=[]):

        # if caltable is a list, insist others are, too
        if (isinstance(caltable,list)):

            if (not isinstance(field,list)):
                field=[field]
            if (not isinstance(intent,list)):
                intent=[intent]
            if (not isinstance(spw,list)):
                spw=[spw]
            if (not isinstance(obs,list)):
                obs=[obs]

            if (not isinstance(tinterp,list)):
                tinterp=[tinterp]
            if (not isinstance(finterp,list)):
                finterp=[finterp]
            if (not isinstance(reach,list)):
                reach=[reach]
            if (not isinstance(calwt,list)):
                calwt=[calwt]

            if (len(obsmap)>0):
                if not isinstance(obsmap[0],list):
                    obsmap=[obsmap]  # nested list
            if (len(fldmap)>0):
                if not isinstance(fldmap[0],list):
                    fldmap=[fldmap]  # nested list
            if (len(spwmap)>0):
                if not isinstance(spwmap[0],list):
                    spwmap=[spwmap]  # nested list
            if (len(antmap)>0):
                if not isinstance(antmap[0],list):
                    antmap=[antmap]  # nested list
            
            igt=0
            for ct in caltable:
                self.parsetorec(caltable=ct,
                                field=field[igt] if (igt<len(field)) else "",
                                intent=intent[igt] if (igt<len(intent)) else "",
                                spw=spw[igt] if (igt<len(spw)) else "",
                                obs=obs[igt] if (igt<len(obs)) else "",
                                tinterp=tinterp[igt] if (igt<len(tinterp)) else "",
                                finterp=finterp[igt] if (igt<len(finterp)) else "",
                                reach=reach[igt] if (igt<len(reach)) else "",
                                calwt=calwt[igt] if (igt<len(calwt)) else calwt[(len(calwt)-1)],
                                obsmap=obsmap[igt] if (igt<len(obsmap)) else [],
                                fldmap=fldmap[igt] if (igt<len(fldmap)) else [],
                                spwmap=spwmap[igt] if (igt<len(spwmap)) else [],
                                antmap=antmap[igt] if (igt<len(antmap)) else [],
                                )
                igt+=1
        else:
            self.parsetorec(caltable=caltable,
                            field=field,intent=intent,spw=spw,obs=obs,
                            tinterp=tinterp,finterp=finterp,
                            reach=reach,calwt=calwt,
                            obsmap=obsmap,fldmap=fldmap,spwmap=spwmap,antmap=antmap)
            
    def parsetorec(self,caltable,
                   field='',intent='',spw='',obs='',
                   tinterp='linear',finterp='',reach='',calwt=True,
                   obsmap=[],fldmap=[],spwmap=[],antmap=[]):

        d0=self.cldefinstance()
        d0["field"]=field
        d0["intent"]=intent
        d0["spw"]=spw
        d0["obs"]=obs

        d0["obsmap"]=obsmap
        d0["fldmap"]=fldmap
        d0["spwmap"]=spwmap
        d0["antmap"]=antmap

        d0["tinterp"]=tinterp
        d0["finterp"]=finterp
        d0["reach"]=reach

        self.addrec({caltable : d0},calwt)

    def addrec(self,crec,calwt):

        ctname=list(crec.keys())[0]

        irec=0
        if (ctname in self.cld):
            # ctname exists, will add a new instance
            irec=len(self.cld[ctname])-1

            # prefer already-set calwt
            calwt0=self.cld[ctname]["calwt"]
            if calwt!=calwt0:
                print('WARNING: For caltable=\''+ctname+'\' using already-set calwt='+str(calwt0)+'.')
        else:
            # ctname does not yet exist, add it
            self.cld[ctname] = {}
            self.cld[ctname]["calwt"]=calwt


        self.cld[ctname][str(irec)]=crec[ctname]

    def list(self):
        print('There are '+str(len(self.cld))+' caltables in the cal library:')
        keys=list(self.cld.keys())
        keys.sort()
        for ct in keys:
            print(ct+': calwt='+str(self.cld[ct]['calwt'])+str(' (')+str(len(self.cld[ct])-1)+str(' instance[s]):'))
            for ims in list(self.cld[ct].keys()):
                if (isinstance(self.cld[ct][ims],dict)):
                    print(' field=\''+str(self.cld[ct][ims]['field'])+'\'', end=' ')
                    print(' intent=\''+str(self.cld[ct][ims]['intent'])+'\'', end=' ')
                    print(' spw=\''+str(self.cld[ct][ims]['spw'])+'\'', end=' ')
                    print(' obs=\''+str(self.cld[ct][ims]['obs'])+'\'')
                    print('  tinterp=\''+str(self.cld[ct][ims]['tinterp'])+'\'', end=' ')
                    print(' finterp=\''+str(self.cld[ct][ims]['finterp'])+'\'')
                    #print ' reach=\''+str(self.cld[ct][ims]['reach'])+'\''
                    print('  obsmap='+str(self.cld[ct][ims]['obsmap']), end=' ')
                    print(' fldmap='+str(self.cld[ct][ims]['fldmap']), end=' ')
                    print(' spwmap='+str(self.cld[ct][ims]['spwmap']), end=' ')
                    print(' antmap='+str(self.cld[ct][ims]['antmap']))


    def write(self,filename,append=False):
        if len(filename)<1:
            raise Exception('Please specify a filename')
        if len(self.cld)<1:
            raise Exception('There is no cal library to write')

        fw="w"
        if append:
            fw="a"
        
        f=open(filename,fw)
        keys0=list(self.cld.keys())
        keys0.sort()
        for ct in keys0:
            ict0=self.cld[ct]
            keys1=list(ict0.keys())
            keys1.sort()
            for ims in keys1:
                ict1=ict0[ims]
                if isinstance(ict1,dict):
                    print('caltable=\''+ct+'\'', end=' ', file=f)
                    print('calwt='+str(ict0['calwt']), end=' ', file=f)
                    if len(ict1['field'])>0:
                        print('field=\''+str(ict1['field'])+'\'', end=' ', file=f)
                    if len(ict1['intent'])>0:
                        print('intent=\''+str(ict1['intent'])+'\'', end=' ', file=f)
                    if len(ict1['spw'])>0:
                        print('spw=\''+str(ict1['spw'])+'\'', end=' ', file=f)
                    if len(ict1['obs'])>0:
                        print('obs=\''+str(ict1['obs'])+'\'', end=' ', file=f)

                    if len(ict1['tinterp'])>0:
                        print('tinterp=\''+str(ict1['tinterp'])+'\'', end=' ', file=f)
                    if len(ict1['finterp'])>0:
                        print('finterp=\''+str(ict1['finterp'])+'\'', end=' ', file=f)
                    if len(ict1['reach'])>0:
                        print('reach=\''+str(ict1['reach'])+'\'', end=' ', file=f)

                    if len(ict1['obsmap'])>0:
                        print('obsmap='+str(ict1['obsmap']), end=' ', file=f)
                    if len(ict1['fldmap'])>0:
                        if isinstance(ict1['fldmap'],str):
                            print('fldmap=\''+str(ict1['fldmap'])+'\'', end=' ', file=f)
                        else:
                            print('fldmap='+str(ict1['fldmap']), end=' ', file=f)
                    if len(ict1['spwmap'])>0:
                        print('spwmap='+str(ict1['spwmap']), end=' ', file=f)
                    if len(ict1['antmap'])>0:
                        print('antmap='+str(ict1['antmap']), end=' ', file=f)

                    print('', file=f)

        f.close()

    def read(self,callibr):

        lines=[]
        if isinstance(callibr,list):
            # a python list of lines has been specified
            lines=callibr
        else:
            # assume a filename has been specified
            lines=open(callibr)

        for line in lines:
            line2=line.strip()  # remove leading/trailing whitespace

            # Attempt to parse if it has content
            if len(line2)>0:
                if line2[0]=='#':
                    # Ignore lines that are comments (or turned off with #)
                    print('Found comment (not parsed): ',line2)
                else:
                    # A nominally parsable line, apparently

                    # reduce whitespace to only singles
                    line2=' '.join(line2.split())  

                    # absorb remaining spaces adjacent to =
                    line2=line2.replace(' =','=')
                    line2=line2.replace('= ','=')

                    # sub , for spaces to delimit keys in the parsed command
                    line2=line2.replace(' ',',')
                    line2=line2.replace(',,',',') # ~corrects likely comma replacement within quotes

                    # add parsetorec() command syntax
                    parsecmd='self.parsetorec('+line2+')'
                    
                    # cope with bool recognition for calwt
                    parsecmd=parsecmd.replace('calwt=T,','calwt=True,') 
                    parsecmd=parsecmd.replace('calwt=T)','calwt=True)')
                    parsecmd=parsecmd.replace('calwt=F,','calwt=False,')
                    parsecmd=parsecmd.replace('calwt=F)','calwt=False)')
                    
                    # execute, and trap/report any errors that occur
                    try:
                        exec(parsecmd)
                    except Exception as errline:
                        self.clear()
                        print('Error: ',errline)
                        print('Problem parsing cal library line (check for typos): "'+line+'"')
                        raise Exception('Problem parsing cal library line (check for typos): '+line)

    def compare(self,other):
        return self.cld==other.cld


def applycaltocallib(filename,append=False,field='',spw='',intent='',
                     gaintable='',gainfield='',interp='',spwmap=[],calwt=True):
    
    if len(filename)<1:
        raise Exception('Please specify a filename')

    if len(gaintable)<1:
        raise Exception('No caltable specified in gaintable')

    c=callibrary()
    c.addold(field=field,spw=spw,intent=intent,gaintable=gaintable,
             gainfield=gainfield,interp=interp,spwmap=spwmap,calwt=calwt)
    c.write(filename,append)
    c.clear()


def testcallib0():

    c=callibrary()
    c.add(caltable="G0",field="0",tinterp="nearest")
    c.add(caltable="G0",field="1,2",tinterp="linear")
    c.add(caltable="B0",tinterp="nearest",finterp="linear")
    c.list()
    return c


def testcallib1():

    c1=callibrary()
    c1.add(caltable=['B','phase','flux'],field='0~1,3~4',
          tinterp=['nearest','linear'],ctfield=['0',''])
    c1.add(caltable=['B','phase','flux'],field='2',
          finterp=['nearest','linear'],ctfield=['0','3,4'])
    print('')
    print('c1: -----------------------------------------------')
    c1.list()

    c2=callibrary()
    c2.add(gaintable=['B','flux'],field='0~4',
          tinterp=['nearest'],ctfield=['0'])
    c2.add(gaintable='phase',field='0~1,3~4',
          tinterp='linear',ctfield='')
    c2.add(gaintable='phase',field='2',
          tinterp='linear',ctfield='3,4')
    print('')
    print('c2: -----------------------------------------------')
    c2.list()

    return (c1,c2)


def testcallib2():

    c1=callibrary()
    c1.add(caltable=['B','flux'],field='0~4',
          tinterp=['nearest'],ctfield=['0'])
    c1.add(caltable='phase',field='0~1,3~4',
          tinterp='linear',ctfield='')
    c1.add(caltable='phase',field='2',
          tinterp='linear',ctfield='3,4')
    c1.write('testcallib2.txt')

    c2=callibrary()
    c2.read('testcallib2.txt')

    print('Cal libraries match?', c2.cld==c1.cld)

    return (c1,c2)


