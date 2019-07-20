#UNFINISHED!!! DO NOT RUN

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
    results=pool.starmap(doBLS,transpose([[time]*nsamp,[mag]*nsamp,[wt]*nsamp,\
                                      makeArrOfCopies(binWt,nsamp),\
                                      makeArrOfCopies(binMag,nsamp),\
                                      period,[nbins]*nsamp,[binExt]*nsamp,\
                                      [minBins]*nsamp,[minWt]*nsamp,\
                                      [totalWt]*nsamp]))
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
    results=pool.starmap(doLS,transpose([[time]*nsamp,[mag]*nsamp,period]))
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

def phaseLightCurve(time,mag,mySmooth,boxSize,p):
    ndata=len(time)
    mySort=[]
    myPhase=[0]*ndata
    myMag=[0]*ndata
    myChi=[0]*ndata
    for j in range(ndata):
        mySort.append([mod(time[j],p)/p,j])
    mySort.sort()
    for j in range(ndata):
        myIdx=mySort[j][1]
        myPhase[j]=mySort[j][0]
        myMag[j]=mag[myIdx]
        if mySmooth:
            mySmooth[j]=-1
        if j>0 and myPhase[j]<myPhase[j-1]:
            print("Sort error?")
    if mySmooth:
        #Now smooth the phases:
        #
        #To reduce computation time, remember values from the
        #last round. Specifically:
        
        #Given that phase[j]>=phase[j-1] (data is sorted)

        #If
        #(phase[j-1] - phase[prevLo-1]) > smooth/2 then
        #(phase[j] - phase[prevLo-1]) > smooth/2 and, generally,
        #(phase[j] - phase[k]) > smooth/2 for any k<prevLo

        #-> there's no need to consider any k<prevLo as the lower
        #edge of the box.

        #If
        #(phase[prevHi] - phase[j-1]) <= smooth/2 then
        #(phase[prevHi] - phase[j]) <= smooth/2 and, generally,
        #(phase[k] - phase[j]) <= smooth/2 for any k <= prevHi

        #-> there's no need to consider any k<prevHi as the upper
        #edge of the box

        ###

        #"Wrapping":
        #Phase-smoothing can still take place for values at the beginning
        #or end of a period by "wrapping" around to the other end of the
        #array. For example, we could smmoth values for phase 0 with
        #those at phase 1-s, or values at phase 1 with those at
        #phase 0+s. Since the assumption in phase-folding is that the
        #signal is periodic, looking at the other end of the array is like
        #wrapping around to the "next" or "previous" period
        bLo=0;bHi=0;prevLo=0;prevHi=0;count=0
        boxSum=0.0;halfbox=boxSize/2.0
        
        for j in range(ndata):
            if PHASE_WRAPPING:
               #Determine value for bHi. If j=0, prevHi=0,
                #thereafter updated to the largest index such
                #that phase[bHi] - phase[j-1]<halfBox
                for bHi in range(prevHi,2*ndata-1):
                    if bHi<(ndata-1):
                        if (myPhase[bHi+1]-myPhase[j])>halfbox:
                            break
                    else:
                        #adjust index and phase range if we're off the
                        #end of the array
                        if (myPhase[bHi-ndata+1]-myPhase[j]+1)>halfbox:
                            break
                if (j==0):
                    #Get initial value for bLo for j=0: step back until
                    #ndata + (bLo -1) is out of the range
                    for bLo in range(0,-ndata,-1): #+1 removed because of >=
                        if (myPhase[j]-myPhase[ndata+bLo-1]+1)>halfbox:
                            break
                else:
                    #once prevLo hase ben properly initialized (i.e. j>0),
                    #start from prevLo and shift box as needed
                    for bLo in range(prevLo,j):
                        if bLo<0:
                            if (myPhase[j]-myPhase[ndata+bLo]+1)<halfbox:
                                break
                        else:
                            if (myPhase[j]-myPhase[bLo])<halfbox:
                                break
            else:
                #Find edges of box (no phase-wrapping)
                for bLo in range(prevLo,j):
                    if (myPhase[j]-myPhase[bLo])<halfbox:
                        break
                for bHi in range(prevHi,ndata-1):
                    if (myPhase[bHi+1]-myPhase[j])>halfbox:
                        break
                if DEBUG:
                    #check that the values for bLo and bHi satisfy
                    #the conditions we were trying to meet
                    if ((myPhase[bHi]-myPhase[j])>halfbox) or\
                       (((bHi+1)<ndata) and ((myPhase[bHi+1]-myPhase[j])<=halfbox))\
                       or ((myPhase[j]-myPhase[bLo])>halfbox) or\
                       ((bLo>0)and((myPhase[j]-myPhase[bLo-1])<=halfbox)):
                        print("BOX ERROR!")

            #Initialize the sum of magnitudes in our box. We will
            #start fro scratch if j = 0 or if our new low edge is
            #above our previous high edge
            if (j==0) or (bLo>=prevHi):
                boxSum=0
                count=0
                for k in range(bLo,bHi+1):#shifted +1 due to <=
                    if k<0:
                        myIdx=ndata+k
                    elif k>=ndata:
                        myIdx=k-ndata
                    else:
                        myIdx=k
                    boxSum+=myMag[myIdx]
                    count+=1
            else:
                #if there is overlap between this box and the
                #previous one, subtract off the left edge and
                #add the right
                for k in range(prevLo,bLo):
                    if k<0:
                        myIdx=ndata+k
                    else:
                        myIdx=k
                    boxSum-=myMag[myIdx]
                    count-=1
                for k in range(prevHi+1,bHi+1):#shifted+= due to <=
                    if k>=ndata:
                        myIdx=k-ndata
                    else:
                        myIdx=k
                    boxSum+=myMag[myIdx]
                    count+=1

            #save the values of bLo and bHi for the next value of j
            prevLo=bLo
            prevHi=bHi
            if count>0:
                mySmooth[j]=boxSum/count
                if myChi:
                    myChi[j]=(myMag[j]-mySmooth[j])**2
    return myChi

def computePlavchan(data,args,fargs,pool):
    #Check for type errors
    if not isinstance(data,dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs,funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args,pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')

    [ndata,time]=funcArgsGetTime(fargs)
    [ndata,mag]=funcArgsGetMag(fargs)
    [data,smooth]=funcArgsGetSmoothedMag(fargs)
    [nsamp,period]=funcArgsGetPeriods(fargs)
    [nsamp,power]=funcArgsGetPower(fargs)
    noutliers=args.nout

    #array to hold the deviation from the smoothed curve for each
    #data point
    [ndata,tmpChi]=funcArgsGetChi(fargs)

    #make sure we don't have more outliers than we have data points
    if noutliers>ndata:
        noutliers=ndata

    #determine reference deviations (recycle "tmpChi" array)
    meanMag=dtGetMean(data,DATA_FIELD_TYPE.DATA_Y)
    for j in range(ndata):
        tmpChi[j]=(mag[j]-meanMag)**2
    tmpChi.sort()

    #sum the values most _poorly_ fit by the model mag=meanMag
    maxStd=0.0;maxChi=0.0 #maxChi is the analogous var for each pd
    for j in range(ndata-1,ndata-noutliers-1,-1):
        maxStd+=tmpChi[j]
    maxStd/=noutliers

    boxSize=fargs.boxSize
    
    #Compute periodogram
    fargs.power=pool.starmap(doPlav,transpose([[time]*nsamp,[mag]*nsamp,\
                                           makeArrOfCopies(smooth,nsamp)\
                                           [boxSize]*nsamp,period]))

def doPlav(time,mag,mySmooth,boxSize,p):
    errval=0
    tmpChi=phaseLightCurve(time,mag,mySmooth,boxSize,maxStd,p)
    tmpChi.sort()
    count=0
    maxChi=0
    for j in range(ndata-1,-1,-1):
        if tmpChi[j]!=errval:
            maxChi+=tmpChi[j]
            count+=1
            if count>=noutliers:
                break
    maxChi/=count
    if maxChi>0:
        return maxStd/maxStd
    return errval

def multiFunc(a,b,c):
    return a/b,c
if __name__=='__main__':
    with mp.Pool(4) as pool:
        results=pool.starmap(multiFunc,transpose([range(10),range(1,11),range(2,12)]))
    
