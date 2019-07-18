#UNFINISHED!!! DO NOT TRY TO RUN!!!
#assumes rectangular matrix
def transpose(l):
    n=[]
    for i in range(len(l)):
        for j in range(len(l[0])):
            if j==len(n):
                n.append([l[i][j]])
            else:
                n[j].append(l[i][j])
    return n
def makeArrOfCopies(inList,num):
    bigArr=[[]]*num
    for k in inList:
        for i in range(len(bigArr)):
            bigArr[i].append(k)
    return bigArr

import multiprocessing as mp
import time
def computeBLS(data,args,fargs,pool):
    #Check for type errors
    if not isinstance(data,dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs,funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args,pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    if not isinstance(pool,mp.Pool):
        raise TypeError('Inappropriate argument type:\npool must be of type Pool!')
    [ndata,time]=funcArgsGetTime(fargs)
    [ndata,mag]=funcArgsGetMag(fargs)
    [nsamp,period]=funcArgsGetPeriods(fargs)
    [nsamp,power]=funcArgsGetPower(fargs)

    #Initializes variables with nonsense values
    blsR=[0.0]*nsamp
    blsS=[0.0]*nsamp
    lowBin0=[0]*nsamp
    lowBin1=[0]*nsamp

    wt=[]
    #In case we want weight as a function of uncertainty
    if WEIGHT_BY_ERR:
        err=dtGetFilteredArray(data,DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY)
    totalWt=0
    for j in range(ndata):
        if WEIGHT_BY_ERR:
            wt.append(err[j])
        else:
            #If we don't, weight everthing the same
            wt.append(1)
        totalWt+=wt[j]

    #If unset, make a value for nbins, based off of ndata
    if args.nbins==DEFAULT_NBINS:
        if ndata<=500:
            args.nbins=50
        elif ndata<=20000:
            args.nbins=int(ndata/10)
        else:
            args.nbins=2000
    nbins=args.nbins

    #Write nbins value to debugfp
    if args.debugfp!=None:
        args.debugfp.write("IN COMPUTE BLS: nbins = "+str(nbins)+"\n")

    #Run checks on FractionOfPeriodInTransitMin/Max before
    #running the algorithm, in case it somehow escaped our
    #pgramArgs.populate() checks
    qmin=args.qmin
    qmax=args.qmax
    if qmin<=0 or qmax<=0 or qmax<qmin:
        raise ValueError("Error: invalid values for qmin/qmax!")

    minBins=int(qmin*nbins)
    if minBins<1:
        minBins=1
        
    minWt=totalWt*qmin
    #I don't love that this is fixed at "5" -- if we convert to
    #weighting with errors, it'll have to change
    if minWt<5:
        minWt=5 #min weight over "low" set of bins

    #maximum number of abins over which a "low" phase can extend:
    #(this is also the amount by which we want to pad the bin array)
    binExt=int(qmax*nbins)+1
    binMax=int(nbins+binExt)

    #Initialize arrays with null values, because we refer to
    #specific elements of these arrays by index later on:
    #just appending as needed won't work here.
    binMag=[0.0]*binMax
    binWt=[0.0]*binMax

    #Compute periodogram
    results=pool.starmap(doBLS,transpose([time]*nsamp,[mag]*nsamp,[wt]*nsamp,\
                                      makeArrOfCopies(binWt,nsamp),\
                                      makeArrOfCopies(binMag,nsamp),\
                                      period,[nbins]*nsamp,[binExt]*nsamp,\
                                      [minBins]*nsamp,[minWt]*nsamp,\
                                      [totalWt]*nsamp))
    results=transpose(results)
    fargs.power=results[0]
    fargs.blsR=results[1]
    fargs.blsS=results[2]
    fargs.lowBin0=results[3]
    fargs.lowBin1=results[4]
#time and mag are arrs. binWt and binMag are working vars,
#so they can't be shared. copies will have to be made for
#each run of the function.
def doBLS(time,mag,wt,binWt,binMag,period,nbins,nsamp,binExt,minBins,minWt,\
          totalWt):
    maxPwr=0.0
    for b in range(nbins):
        binMag[b]=0.0
        binWt[b]=0.0
    ndata=len(time)
    binMax=len(binWt)
    for j in range(ndata):
        phase=(time[j]/period)%1
        b=int(math.floor(nbins*phase))
        binWt[b]+=wt[j]
        binMag[b]+=wt[j]*mag[j]
    for b in range(nbins,binMax):
        binWt[b]=binWt[b-nbins]
        binMag[b]=binMag[b-nbins]
    maxPwr=0.0
    for b in range(nbins):
        binCt=0
        sumWt=0.0
        sumMag=0.0
        for k in range(b,b+binExt+1):
            binCt+=1
            sumWt+=binWt[k]
            sumMag+=binMag[k]
            if binCt>=minBins and sumWt>=minWt and\
               sumWt<totalWt:
                pwr=(sumMag**2)/(sumWt*(totalWt-sumWt))
                if pwr>=maxPwr:
                    maxPwr=pwr
                    lowStart=b
                    lowEnd=k
                    lowWt=sumWt
                    lowMag=sumMag
    maxPwr=math.sqrt(maxPwr)
    if maxPwr>0:
        return (maxPwr,lowWt/totalWt,lowMag,lowStart,lowEnd)
    return (0,0,0,0,0)
def computeLombScargle(data,args,fargs,pool):
    if not isinstance(pool,mp.Pool):
        raise TypeError('Inappropriate argument type:\npool must be of type Pool!')
    if not isinstance(data,dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs,funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args,pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    [ndata,time]=funcArgsGetTime(fargs)
    [ndata,mag]=funcArgsGetMag(fargs)
    [nsamp,period]=funcArgsGetPeriods(fargs)
    [nsamp,power]=funcArgsGetPower(fargs)

    #Compute stats on magnitude
    sdMag=dtGetDev(data,DATA_FIELD_TYPE.DATA_Y)
    if sdMag==0:
        raise ValueError("Error in InputFile: Zero deviation in data values!")
    results=pool.starmap(doLS,transpose([time]*ndata,[mag]*ndata,period))
    fargs.power=results

def doLS(time,mag,p):
    ndata=len(time)
    w=2*math.pi/p
    tnum=0
    tdenom=0
    for j in range(ndata):
        tnum+=math.sin(2.0*w*time[j])
        tdenom+=math.cos(2.0*w*time[j])
    t=(1/(2*w))*math.atan2(tnum,tdenom)
    lnum=0
    ldenom=0
    rnum=0
    rdenom=0
    for j in range(ndata):
        s=math.sin(w*(time[j]-t))
        c=math.cos(w*(time[j]-t))
        rnum+=mag[j]*s
        lnum+=mag[j]*c
        rdenom+=s*s
        ldenom+=c*c
    return (1/(2*sdMag*sdMag))*((lnum*lnum)/ldenom+(rnum*rnum/rdenom))

def multiFunc(a,b,c):
    return a/b,c
if __name__=='__main__':
    with mp.Pool(4) as pool:
        results=pool.starmap(multiFunc,transpose([range(10),range(1,11),range(2,12)]))
    

