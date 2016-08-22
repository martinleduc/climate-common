import os        
def simget(machines, sims,datadir='/home/leduc/data/CCR/',gettsi=1):
    '''Importer localement les series temporelles de simulation generees sur differentes machines.'''
    i=0
    for sim in sims.split():
        machine=machines.split()[i]
        print ''
        print '.......... Import '+sim+' from '+machine+' ..........'

        if machine=='climate1':
            remotedir='models/RUN/'
            remoteloc="-e ssh m_ledu@climate-1.concordia.ca:"+remotedir+sim
        elif machine=='pharos1':
            remotedir='/exec1/leduc/'
            remoteloc="pharos1:"+remotedir+sim

        elif machine=='pharos2':
            remotedir='/exec2/leduc/'
            remoteloc="pharos2:"+remotedir+sim
                        
        elif machine=='pharos3':
            remotedir='/exec3/leduc/RUN/'
            remoteloc="-e ssh pharos3:"+remotedir+sim
        elif machine=='ubuntu':
            remotedir='/home/leduc/UVic_ESCM/RUN/'
            remoteloc=remotedir+sim
            
        if gettsi==1:
            if machine[0:6]=='pharos':
                os.system("rcp "+remoteloc+"/tsi*.nc"+" "+datadir+"/"+sim+"/")                
            else:
                os.system("rsync -vzrlptD "+remoteloc+"/tsi*.nc"+" "+datadir+"/"+sim+"/")
        elif gettsi==0:
            if machine[0:6]=='pharos':
                os.system("rcp -r "+remoteloc+"/*"+" "+datadir+"/"+sim+"/")                
            else:            
                os.system("rsync -vzrlptD "+remoteloc+" "+datadir)

        i+=1
    #pm.plottsi('A_sat A_totcarb L_totcarb O_totcarb',sims,datadir,outdir=outdir)
         
        
    
