#!/usr/bin/env python
import sys
import numpy
import array

import pywcs

import math
import pyfits
import os

import ROOT


#check for matplotlib
#HASPYLAB = True
#try:
#    import pylab, matplotlib
#    #rc('backend',TkAgg)
#    try:
#        dobatch = os.getenv("PYLAB_BATCH")
#        if dobatch == 'True':
#            matplotlib.interactive(False)
#            print 'PYLAB WILL RUN IN BATCH MODE'
#    except NameError:
#        print 'PYLAB WILL RUN INTERACTIVELY'
#except:
HASPYLAB = False
#print 'HASPYLAB= ',HASPYLAB
##def plotLC(det):
##    if isintance(det,GBM.GBM):
##        pass
##    else:
##        pass
##    f=pyfits.open('/home/cohen/data1/WORK/GRBSCRIPTS/data/OUT/090103597/090103597_N1_lc.fits')
##    t=f[1].data.field('TIME')-f[1].header['TSTART']
##    c=f[1].data.field('COUNTS')
##    #matplotlib!  plot(t,c)



def PlotTS(tsmap,quick=False,markers=[]):
    '''
    Code by Frederic Piron (piron@in2p3.fr):
    - reads & plots the TS maps produced by gttsmap
    - computes equivalent error radii at different confidence levels 
    
    pywcs code part adapted from "Loading WCS information from a FITS file":
    http://stsdas.stsci.edu/astrolib/pywcs/examples.html#loading-wcs-information-from-a-fits-file
    pywcs doc can be found here:
    http://stsdas.stsci.edu/astrolib/pywcs/
    '''
    
    ############################################
    ### initializations
    ############################################
    # prepare palette for B&W plot
    nRGBs = 2
    ncont = 100
    stops = [ 0.00,  1.00 ]#100% from dark grey to white
    red   = [ 0.30,  1.00 ]
    green = [ 0.30,  1.00 ]
    blue  = [ 0.30,  1.00 ]
    stopsArray = array.array('d', stops)
    redArray   = array.array('d', red)
    greenArray = array.array('d', green)
    blueArray  = array.array('d', blue)
    mycol=ROOT.TColor
    FI=mycol.CreateGradientColorTable(nRGBs, stopsArray, redArray, greenArray, blueArray, ncont)
    colors=[]
    for i in range(ncont):
        colors.append(FI+i)
        colorArray = array.array('i', colors)
        pass
    #ROOT.gStyle.SetNumberContours(50)#uncomment this line to make image contours very smooth
    
    # ROOT style
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)

    # pi/180
    pio180=math.pi/180.

    # create file name for final plot
    fout=tsmap.replace('.fits','.png')
    foutbw=tsmap.replace('.fits','BW.png')
    fout_eps=tsmap.replace('.fits','.eps')
    foutbw_eps=tsmap.replace('.fits','BW.eps')
    ############################################
    ### load the FITS hdulist, parse the WCS keywords in the primary HDU and get data
    ############################################
    print '<><><><> reading the TS map...'
    hdulist = pyfits.open(tsmap)
    mywcs=pywcs.WCS(hdulist[0].header,key=' ')
    data=hdulist[0].data

    ############################################
    ### store TS map in 2D histogram with same binning / center point
    ############################################
    nbx=mywcs.naxis1
    nby=mywcs.naxis2
    deltax=mywcs.wcs.cdelt[0]
    deltay=mywcs.wcs.cdelt[1]
    xref=mywcs.wcs.crpix[0]
    yref=mywcs.wcs.crpix[1]
    raref=mywcs.wcs.crval[0]
    decref=mywcs.wcs.crval[1]
    print '%d x %d pixellised map with %6.3f deg x %6.3f deg binning, centered at (x, y)=(%6.3f, %6.3f) or (RA, Dec)=(%7.3f, %7.3f)' %(nbx,nby,deltax,deltay,xref,yref,raref,decref)
    xmin=xref-nbx/2.
    xmax=xref+nbx/2.
    ymin=yref-nby/2.
    ymax=yref+nby/2.
    hTS=ROOT.TH2F('hTS','',nbx,xmin,xmax,nby,ymin,ymax)
    for i in range(nbx):
        for j in range(nby):
            hTS.SetBinContent(i+1,j+1,data[j][i])#bin 0 of histogram is for underflow
            pass
        pass
    ROOT.SetOwnership(hTS,False)

    canTSbw=ROOT.TCanvas('canTSbw','',0,0,600,600)
    canTS=ROOT.TCanvas('canTS','',0,600,600,600)
    ROOT.SetOwnership(canTS,False)
    ROOT.SetOwnership(canTSbw,False)
    
    canTS.cd()
    hTS.Draw('ah*colz')
    hTS.SetLabelFont(42,"xyz");
    hTS.SetLabelSize(0.02,"xyz");
    canTS.Update()

    ############################################
    ### find maximum
    ############################################
    grbi=ROOT.Long(0)
    grbj=ROOT.Long(0)
    kkk=ROOT.Long(0)
    hTS.GetMaximumBin(grbi,grbj,kkk)
    grbts=hTS.GetBinContent(grbi,grbj)
    grbx=hTS.GetXaxis().GetBinCenter(grbi)
    grby=hTS.GetYaxis().GetBinCenter(grbj)
    grbpix=numpy.array([[grbx,grby]],numpy.float_)
    grbsky=mywcs.all_pix2sky(grbpix,1)
    grbra=grbsky[0,0]
    grbdec=grbsky[0,1]
    print "Maximum TS=%7.3f in bin (%d, %d) at (x, y)=(%6.3f, %6.3f) or (RA, Dec)=(%7.3f, %7.3f)" %(grbts,grbi,grbj,grbx,grby,grbra,grbdec)

    ############################################
    ### create map of deltaTS and define its contours
    ############################################
    hDTS=ROOT.TH2F('hDTS','',nbx,xmin,xmax,nby,ymin,ymax)
    for i in range(nbx):
        for j in range(nby):
            hDTS.SetBinContent(i+1,j+1,grbts-data[j][i])
            pass
        pass
    ROOT.SetOwnership(hDTS,False)
    if quick:
        deltaTS=[ 2.30, 4.61]#delta_TS values (sorted in increasing order)
        CLprob =[   68,   90]#C.L. levels to probe
        CLdisp =[ True, True]#if false, compute error but do not display
        pass
    else:
        deltaTS=[ 2.30, 4.61,  5.99,  9.21]#delta_TS values (sorted in increasing order)
        CLprob =[   68,   90,    95,    99]#C.L. levels to probe
        CLdisp =[ True, True, False, False]#if false, compute error but do not display
        pass

    hDTS.SetContour(len(CLprob))
    for il in range(len(CLprob)):
        hDTS.SetContourLevel(il,deltaTS[il])
        pass

    ############################################
    ### retrieve the contour graphs and superimpose them on TS map with cross at maximum
    ############################################
    print '<><><><> retrieving and plotting the contour graphs...'
    canTS.cd()
    hDTS.Draw('cont list')#the list options create the list of contour graphs
    canTS.Update()
    hTS.Draw('ah*colz')
    icol=1
    grbpos=ROOT.TMarker(grbx,grby,3)
    grbpos.SetMarkerSize(1)
    grbpos.SetMarkerColor(icol)
    grbpos.Draw()    
    ROOT.SetOwnership(grbpos,False)
    canTS.Update()
    contours=ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    tgcont=[]#tgraphs of C.L. contour in TS map coordinates
    CLprob2=[]#final C.L. levels (some contours might not be found)
    CLdisp2=[]
    area=[]#area of the C.L. contour in RA, Dec coordinates; will be used to compute the final error radii
    #print "Total number of contours: %d" %(contours.GetSize())
    #for il in range(contours.GetSize()):
        #LevelContList=contours.At(il)
        #LevelContNb=LevelContList.GetSize()
        #print "Contour %d has %d graph(s)" %(il,LevelContNb)

    for il in range(contours.GetSize()):
        # add an entry in tables
        area.append(0.)
        CLprob2.append(CLprob[il])
        CLdisp2.append(CLdisp[il])
        # retrieve contour list for the current confidence level
        LevelContList=contours.At(il)
        LevelContNb=LevelContList.GetSize()
        print "<><> %d percent level contour has %d graph(s)" %(CLprob[il],LevelContNb)
        # retrieve first contour
        tgcont.append(ROOT.TGraph())
        tgcont[-1]=LevelContList.First()
        # for this confidence level, look for the contour which contains the GRB position
        foundgraph=False
        for ig in range (LevelContNb):
            if (ig!=0):
                tgcont[-1]=LevelContList.After(tgcont[-1])
            np=tgcont[-1].GetN()
            xp=[]
            yp=[]
            xp=tgcont[-1].GetX()
            yp=tgcont[-1].GetY()
            print "Graph %d has %d points starting at (x, y)= (%6.3f, %6.3f)" %(ig,np,xp[0],yp[0])
            if ((xp[0]!=xp[np-1]) or (yp[0]!=yp[np-1])):
                print "%d percent contour is not closed (contained in the map): error won't be accurate" %(CLprob[il])
            j=np-1
            oddNodes=False
            for i in range(np):
                if ((yp[i]<grby and yp[j]>=grby) or (yp[j]<grby and yp[i]>=grby)):
                    if (xp[i]+(grby-yp[i])/(yp[j]-yp[i])*(xp[j]-xp[i])<grbx):
                        oddNodes=not oddNodes
                        pass
                    pass
                j=i
            if (not oddNodes):#the position does not fall within the contour
                print "Graph %d does not contain the GRB position -> ignored" %(ig)
                pass
            else:
                print "Graph %d contains the GRB position -> let's use that one (and ignore other graphs if any)" %(ig)
                foundgraph=True
                break
        # found graph: set style and draw (if in the list of contours to draw)
        if (foundgraph):
            tgcont[-1].SetLineWidth(1)
            tgcont[-1].SetLineStyle(il+1)
            if CLdisp2[-1]:
                tgcont[-1].Draw("sameL")
                ROOT.SetOwnership(tgcont[-1],False)
            pass
        else:
            print '%d percent C.L. contour was not found' %(CLprob[il])
            tgcont.pop()
            area.pop()
            CLprob2.pop()
            CLdisp2.pop()
            pass
        pass
    canTS.Update()


    # idem in B&W
    canTSbw.cd()
    hTS.Draw('ah*colz')
    grbpos.Draw()    
    for il in range(len(tgcont)):
        if CLdisp2[il]:
            tgcont[il].Draw()
        pass
    canTSbw.Update()

    if len(tgcont)==0:
        print '******** no contour was found: exiting'
        canTS.Print(fout)
        canTS.Print(fout_eps)
        ROOT.gStyle.SetPalette(ncont,colorArray)
        canTSbw.Print(foutbw)
        canTSbw.Print(foutbw_eps)
        ROOT.gStyle.SetPalette(1)
        return grbts,grbra,grbdec,0.,0.,0.,0.
        pass

    ############################################
    ### get boundaries of (RA,Dec) map
    ############################################
    print '<><><><> scanning the TS map to get RA and Dec boundaries...'
    raArr=[]
    decArr=[]
    if quick:
        nsubx=nbx
        nsuby=nby
        pass
    else:
        nsubx=10*nbx
        nsuby=10*nby
        pass
    for i in range(nsubx+1):
        pixx=xmin+i*(xmax-xmin)/nsubx
        for j in range(nsuby+1):
            pixy=ymin+j*(ymax-ymin)/nsuby
            pix=numpy.array([[pixx,pixy]],numpy.float_)
            sky=mywcs.all_pix2sky(pix,1)
            raArr.append(sky[0,0])
            decArr.append(sky[0,1])
            pass
        pass
    ramin=min(raArr)
    ramax=max(raArr)
    decmin=min(decArr)
    decmax=max(decArr)
    print 'TS map: RA min / max = %7.3f / %7.3f ; Dec min / max = %7.3f / %7.3f' %(ramin,ramax,decmin,decmax)

    ############################################
    ### draw (RA, Dec) lines
    ############################################
    print '<><><><> drawing the (RA, Dec) lines...'

    frac=0.2#(RA, Dec) grid with ~1/frac lines in each direction

    exp=int(math.log10(ramax-ramin))
    scale=math.pow(10.,1.-exp)
    deltara=int(frac*math.ceil(ramax-ramin)*scale)/scale
    
    exp=int(math.log10(decmax-decmin))
    scale=math.pow(10.,1.-exp)
    deltadec=int(frac*math.ceil(decmax-decmin)*scale)/scale

    print 'ramax-ramin=%f decmax-decmin=%f -> deltara=%f deltadec=%f' %(ramax-ramin,decmax-decmin,deltara,deltadec)

    tgline=[]#RA and Dec lines
    tglinelab=[]#labels for these lines

    fracin=0.05
    xminin=xmin+fracin*(xmax-xmin)
    xmaxin=xmax-fracin*(xmax-xmin)
    yminin=ymin+fracin*(ymax-ymin)
    ymaxin=ymax-fracin*(ymax-ymin)

    dec=decref-int((decref-decmin)/deltadec)*deltadec
    while dec<=decmax:
        tgline.append(ROOT.TGraph())
        num=0
        if deltax>0:
            ra=ramin
        else:
            ra=ramax
        nolabelyet=True
        while (deltax>0 and ra<=ramax) or (deltax<0 and ra>=ramin):
            sky=numpy.array([[ra,dec]],numpy.float_)
            pix=mywcs.wcs_sky2pix(sky,1)
            pixx=ROOT.Double(pix[0,0])
            pixy=ROOT.Double(pix[0,1])
            #print "(ra, dec)=(%f, %f) -> (x, y)=(%f, %f)" %(ra,dec,pixx,pixy)
            if pixx>xmin and pixx<xmax and pixy>ymin and pixy<ymax:#does the current label position fall in the map?
                tgline[-1].SetPoint(num,pixx,pixy)
                num+=1
                if nolabelyet and pixx>xminin and pixx<xmaxin and pixy>yminin and pixy<ymaxin:#no label yet and well in the map?
                    nolabelyet=False
                    tglinelab.append(ROOT.TText(pixx,pixy,'%6.2f'%dec))#Dec line label
                    pass
                pass
            if deltax>0:
                ra+=(ramax-ramin)/1000
            else:
                ra-=(ramax-ramin)/1000
                pass
            pass
        if num==0:
            tgline.pop()
            pass
        dec+=deltadec
        pass

    ra=raref-int((raref-ramin)/deltara)*deltara
    while ra<=ramax:
        tgline.append(ROOT.TGraph())
        num=0
        dec=decmin
        nolabelyet=True
        while dec<=decmax:
            sky=numpy.array([[ra,dec]],numpy.float_)
            pix=mywcs.wcs_sky2pix(sky,1)
            pixx=ROOT.Double(pix[0,0])
            pixy=ROOT.Double(pix[0,1])
            if pixx>xmin and pixx<xmax and pixy>ymin and pixy<ymax:#does the current label position fall in the map?
                tgline[-1].SetPoint(num,pixx,pixy)
                num+=1
                if nolabelyet and pixx>xminin and pixx<xmaxin and pixy>yminin and pixy<ymaxin:#no label yet and well in the map?
                    nolabelyet=False
                    tglinelab.append(ROOT.TText(pixx,pixy,'%6.2f'%ra))#RA line label
                    pass
                pass
            dec+=(decmax-decmin)/1000
            pass
        if num==0:
            tgline.pop()
            pass
        ra+=deltara
        pass

    # draw lines
    for i in range(len(tgline)):
        tgline[i].SetLineStyle(3)
        tgline[i].SetLineColor(19)
        ROOT.SetOwnership(tgline[i],False)

        canTS.cd()
        tgline[i].Draw("sameC" )
        canTS.Update()

        canTSbw.cd()
        tgline[i].Draw("sameC")
        canTSbw.Update()
        pass

    # draw labels
    for i in range(len(tglinelab)):
        tglinelab[i].SetTextColor(19)
        tglinelab[i].SetTextSize(0.02)
        ROOT.SetOwnership(tglinelab[i],False)
        
        canTS.cd()
        tglinelab[i].Draw()
        canTS.Update()

        canTSbw.cd()
        tglinelab[i].Draw("same")
        canTSbw.Update()
        pass

    ############################################
    ### compute error radii
    ############################################
    print '<><><><> computing error radii...'
    np=[]
    xp=[]
    yp=[]
    errok=[]
    for il in range(len(tgcont)):
        np.append(tgcont[il].GetN())
        xp.append(tgcont[il].GetX())
        yp.append(tgcont[il].GetY())
        if ((xp[il][0]!=xp[il][np[il]-1]) or (yp[il][0]!=yp[il][np[il]-1])):
            print "******** %d percent contour is not closed (contained in the map): error won't be accurate" %(CLprob2[il])
            errok.append(0.)
            pass
        else:
            errok.append(1.)
            pass
        pass
    # RA, Dec step is 0.5% of their range; take 1/5th of pixel size if range is too big (e.g. in RA close to the pole)
    if quick:#quick mode, assuming x=RA and Dec=y, which is ~correct only for small scales and far from the pole
        deltara=math.fabs(deltax)
        deltadec=math.fabs(deltay)
        pass
    else:
        deltara=min(0.005*(ramax-ramin),0.2*math.fabs(deltax))
        deltadec=min(0.005*(decmax-decmin),0.2*math.fabs(deltay))
        pass
    print 'Using (RA, Dec) binning of (%f, %f)' %(deltara,deltadec)

    grin=ROOT.TGraph()#plot this tgraph further below to check if areas are correctly computed, i.e. with a binning which is smaller than the size of a pixel in the map
    nin=0
    ntotiter=int((decmax-decmin)/deltadec)
    dniter=int(0.1*ntotiter)
    niter=0
    # loop over RA,Dec infinitesimal bins
    dec=decmin
    while dec<=decmax:
        niter+=1
        if niter%dniter==0:
            per=niter*100./(ntotiter*1.)
            print '%4.0f percent done...' %(per)
            pass
        ra=ramin
        while ra<=ramax:
            #print "(RA, Dec)=(%7.3f, %7.3f)" %(ra,dec)        
            # compute area for the current RA,Dec infinitesimal bin
            pixarea=deltara*math.cos(dec*pio180)*deltadec
            # get projected position in map for the current RA,Dec infinitesimal bin
            sky=numpy.array([[ra,dec]],numpy.float_)
            pix=mywcs.wcs_sky2pix(sky,1)
            pixx=ROOT.Double(pix[0,0])
            pixy=ROOT.Double(pix[0,1])
            # check if this position lies within the C.L. contours
            for il in range(len(tgcont)):
                #see http://root.cern.ch/root/htmldoc/src/TMath.h.html#vuXOW
                #// Function which returns kTRUE if point x,y lies inside the
                #// polygon defined by the np points in arrays xp and yp
                #// Note that the polygon may be open or closed.
                j=np[il]-1
                oddNodes=False
                for i in range(np[il]):
                    if ((yp[il][i]<pixy and yp[il][j]>=pixy) or (yp[il][j]<pixy and yp[il][i]>=pixy)):
                        if (xp[il][i]+(pixy-yp[il][i])/(yp[il][j]-yp[il][i])*(xp[il][j]-xp[il][i])<pixx):
                            oddNodes=not oddNodes
                            pass
                        pass
                    j=i
                if (oddNodes):#the position falls within the contour
                    #print "(RA, Dec)=(%7.3f, %7.3f) -> (x,y)=(%6.3f, %6.3f) is inside (pixarea=%f)" %(ra,dec,pixx,pixy,pixarea)
                    if (il==0):
                        grin.SetPoint(nin,pixx,pixy)
                        nin+=1
                        pass
                    # add the area of the current RA,Dec infinitesimal bin to the total area
                    area[il]+=pixarea
                    pass
                pass
            ra+=deltara
            pass
        dec+=deltadec
        pass

    grin.SetMarkerStyle(19)
    grin.SetMarkerColor(7)
    grin.SetLineColor(7)

    #canTS.cd()
    #grin.Draw("sameP")#uncomment this line to check area computation
    #ROOT.SetOwnership(grin,False)
    #canTS.Update()

    # idem in B&W
    #canTSbw.cd()
    #grin.Draw("sameP")#uncomment this line to check area computation
    #canTSbw.Update()

    err=[]
    for il in range(len(tgcont)):
        area[il]=area[il]*pio180**2.
        # final error radius is the 1/2 angle of the spherical cap with same area as the C.L. contour
        err.append(math.acos(1.-area[il]/(2.*math.pi))*180./math.pi)
        if errok[il]==1:
            print '%d percent C.L. error = %8.4f deg (area = %8.2e sr)' %(CLprob2[il],err[il],area[il])
        else:
            print '%d percent C.L. error = %8.4f deg **** not accurate **** (area = %8.2e sr)' %(CLprob2[il],err[il],area[il])
            pass
        pass

    ############################################
    ### draw error circles
    ############################################
    print '<><><><> drawing the error circles...'
    # burst direction
    cp = math.cos(grbra*pio180)
    sp = math.sin(grbra*pio180)
    ct = math.cos(grbdec*pio180)
    st = math.sin(grbdec*pio180)
    kx=ct*cp
    ky=ct*sp
    kz=st
    
    tgCcont=[]
    for il in range(len(tgcont)):
        tgCcont.append(ROOT.TGraph())
        n2=0
        # start with direction on same meridian, making an angle err[il] with dir
        if grbdec>0:
            ct = math.cos((grbdec-err[il])*pio180)
            st = math.sin((grbdec-err[il])*pio180)
            pass
        else:
            ct = math.cos((grbdec+err[il])*pio180)
            st = math.sin((grbdec+err[il])*pio180)
            pass
        vx=ct*cp
        vy=ct*sp
        vz=st
        # compute scalar product (should be equal to err[il])
        scalprod=vx*kx+vy*ky+vz*kz
        #print 'check err=%f deg scalprod=%f deg' %(err[il],math.acos(scalprod)/pio180)
        vecprodx = ky*vz-kz*vy
        vecprody = kz*vx-kx*vz
        vecprodz = kx*vy-ky*vx
        # create rotated directions around burst direction
        na=1000
        for ia in range(na+1):
            rotang=2.*math.pi*ia/na 
            crot=math.cos(rotang)
            srot=math.sin(rotang)
            # see Rodrigues' rotation formula
            # http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
            rvx = vx*crot + vecprodx*srot + kx*scalprod*(1.-crot)
            rvy = vy*crot + vecprody*srot + ky*scalprod*(1.-crot)
            rvz = vz*crot + vecprodz*srot + kz*scalprod*(1.-crot)
            if rvx!=0:
                rra = math.atan(rvy/rvx)/pio180
                if rvx<0:
                    rra+=180
                    pass
                if rra<0:
                    rra+=360#RA is between 0 and 360
                    pass
                pass
            else:
                rra=90
                if rvy<0:
                    rra=270
                    pass
                pass
            rdec=math.asin(rvz)/pio180
            # go back to pixel coordinates
            rsky=numpy.array([[rra,rdec]],numpy.float_)
            rpix=mywcs.wcs_sky2pix(rsky,1)
            rpixx=ROOT.Double(rpix[0,0])
            rpixy=ROOT.Double(rpix[0,1])
            #print "rotang=%f (grbra, grbdec)=(%f, %f) (ra, dec)=(%f, %f) -> (x, y)=(%f, %f)" %(rotang,grbra,grbdec,rra,rdec,rpixx,rpixy)
            tgCcont[il].SetPoint(n2,rpixx,rpixy)
            n2+=1
            pass
        tgCcont[il].SetLineColor(icol)
        tgCcont[il].SetLineWidth(2)
        tgCcont[il].SetLineStyle(il+1)
        canTS.cd()
        if CLdisp2[il]:
            tgCcont[il].Draw("sameL")
            ROOT.SetOwnership(tgCcont[il],False)
        canTS.Update()

        # idem in B&W
        canTSbw.cd()
        if CLdisp2[il]:
            tgCcont[il].Draw("sameL")
        canTSbw.Update()
        pass
    
    ############################################
    ### caption
    ############################################
    caption=ROOT.TPaveText(0.05,0.88,0.55,.98,"NDC");
    caption.SetFillColor(0);
    caption.SetBorderSize(1);
    if err[-1]>1.:
        title='Maximum TS=%5.1f at (RA, Dec)=(%6.2f, %6.2f)' %(grbts,grbra,grbdec)
    else:
        title='Maximum TS=%5.1f at (RA, Dec)=(%7.3f, %7.3f)' %(grbts,grbra,grbdec)
        pass
    caption.AddText(title);
    title='Error radius (deg)'
    first=True
    for il in range(len(tgcont)):
        if CLdisp2[il]:
            if first:
                first=False
                title+=': '
                pass
            else:
                title+=', '
                pass
            if err[-1]>1.:
                title+='%5.2f' %(err[il])
            else:
                title+='%6.3f' %(err[il])
                pass
            
            title+=' ('+str(CLprob2[il])+'%'+')'
            pass

    caption.AddText(title);    

    canTS.cd()
    caption.Draw();
    ROOT.SetOwnership(caption,False)
    canTS.Update()
    
    canTSbw.cd()
    caption.Draw()
    canTSbw.Update()

    ############################################
    ### save final plot
    ############################################
    print '<><><><> saving the image in color and in B&W...'
    canTS.Print(fout)
    canTS.Print(fout_eps)
    ROOT.gStyle.SetPalette(ncont,colorArray)
    canTSbw.Print(foutbw)
    canTSbw.Print(foutbw_eps)
    ROOT.gStyle.SetPalette(1)

    ############################################
    ### returned values: best position and errors
    ############################################
    try:
        err68=err[CLprob2.index(68)]
    except:
        err68=0.
    try:
        err90=err[CLprob2.index(90)]
    except:
        err90=0.
    try:
        err95=err[CLprob2.index(95)]
    except:
        err95=0.
    try:
        err99=err[CLprob2.index(99)]
    except:
        err99=0.

    print '<><><><> done!'
    
    return grbts,grbra,grbdec,err68,err90,err95,err99


if __name__=='__main__':
    import sys
    PlotTS(sys.argv[1])
