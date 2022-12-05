'''
# This code compares equivalent cases (command lines)
# for the scattering of Bessel beams calculated in ADDA.
'''

import os,shutil,re,math
import numpy as np
import matplotlib.pyplot as plt



fdiff = 1e-3 # fixed relative difference

pi = 3.1415926535897932384626433832795029


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda")
dirname = 'out'    

                                                                    

def extractCrosSeq(mode,pol):
    dname = dirname + str(mode)
    files = os.listdir(dname)
    with open(dname + '/' + str(files[files.index('CrossSec-' + pol)])) as f:
        file = f.read()
    f.close()
    dat = re.split('[\n \t]',file)
    Cext = float(dat[2])
# =============================================================================
#     Qext = float(dat[5])
#     Cabs = float(dat[8])
#     Qabs = float(dat[11])
# =============================================================================
    
    #return Cext,Qext,Cabs,Qabs
    return Cext

def extractBeam(mode,pol):
    dname = dirname + str(mode)
    files = os.listdir(dname)
    with open(dname + '/' + str(files[files.index('IncBeam-' + pol)])) as f:
        file = f.read()
    f.close()
    dat = re.split('[\n \t]',file)
    E2 = np.array(dat[13:-1:10],dtype='float')
    Exr = np.array(dat[14:-1:10],dtype='float')
    Exi = np.array(dat[15:-1:10],dtype='float')
    Eyr = np.array(dat[16:-1:10],dtype='float')
    Eyi = np.array(dat[17:-1:10],dtype='float')
    Ezr = np.array(dat[18:-1:10],dtype='float')
    Ezi = np.array(dat[19:-1:10],dtype='float')
    
    return E2,Exr,Exi,Eyr,Eyi,Ezr,Ezi

# data generation (run of ADDA code)
def adda_run(mode,option):                                                            
    dname = dirname + str(mode)
    #os.makedirs(dname, exist_ok=True)
    #cmdline = adda_exec + ' -store_beam -dir ' + dname + option
    cmdline = adda_exec + ' -store_beam -eps 8 -dir ' + dname + option + ' > ' + os.devnull
    #cmdline = adda_exec + ' -dir ' + dname + option
    os.system(cmdline)
    
# plot gata on angle
def plotAngle(l1,l2):
    nA = 10
    arrA = np.linspace(10,100,nA)
    E1 = np.zeros(nA)
    E2 = np.zeros(nA)
    CS = np.zeros(nA)
    
    flag = 1
    
    for i in range(nA):
        shutil.rmtree('out1',ignore_errors=True)
        shutil.rmtree('out2',ignore_errors=True)
        #adda_run(1,l1+' '+str(arrA[i]))
        adda_run(1,l1)
        adda_run(2,l2+' '+str(arrA[i]))
        e12,e1xr,e1xi,e1yr,e1yi,e1zr,e1zi = extractBeam(1,'Y')
        e22,e2xr,e2xi,e2yr,e2yi,e2zr,e2zi = extractBeam(2,'Y')
        c1 = extractCrosSeq(1,'Y')
        c2 = extractCrosSeq(2,'Y')
        
        if flag == 0:
            E1[i] = np.mean(e12)
            E2[i] = np.mean(e22)
        #else: E1[i] = np.mean(np.abs(e12-e22)/e12)
        else: E1[i] = np.max(np.abs(((e1xr-e2xr+e1yr-e2yr+e1zr-e2zr)+1j*(e1xi-e2xi+e1yi-e2yi+e1zi-e2zi))/e12))
        #else: E1[i] = np.mean(np.abs((e1xr-e2xr+e1yr-e2yr+e1zr-e2zr)+1j*(e1xi-e2xi+e1yi-e2yi+e1zi-e2zi)))
        #else: E1[i] = np.mean(np.abs(((e1zr-e2zr)+1j*(e1zi-e2zi))/(e1zr+1j*e1zi)))
        #else: E1[i] = np.mean(np.abs(((e1zi-e2zi))/(e2zi)))
        #else: E1[i] = np.max(np.abs((e1zr-e2zr)+1j*(e1zi-e2zi))/e22)
        
        CS[i] = abs((c1-c2)/c1)
        
    fig = plt.figure(figsize=(9.5,5))
    ax = fig.add_subplot(121)
    #rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    if flag == 0:
        plt.plot(arrA, E1, label = 'line 1', color = 'b')
        plt.plot(arrA, E2, label = 'line 2', color = 'm',linestyle='dashed')
        plt.legend()
        plt.ylabel(r'$I$')
    else:
        plt.plot(arrA, E1, label = 'line 1', color = 'r', marker='.', linewidth=0.5)
        plt.ylabel(r'$I$ max relative error')
    plt.yscale('log')
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Half-cone angle $\alpha$, deg')
    
    #fig = plt.figure(figsize=(6.5,5))
    ax = fig.add_subplot(122)
    lev = fdiff*np.ones(nA)
    plt.plot(arrA, CS, color = 'k', marker='.', linewidth=0.5)
    plt.plot(arrA, lev, color = 'r', linewidth=0.5)
    plt.yscale('log')
    plt.ylabel(r'$Cext$ relative error')
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Half-cone angle $\alpha$, deg')

    plt.tight_layout()
    plt.show()
    
        

        
comm = ' -grid 16 '
line1 = comm + ' -beam davis3   10 '
line2 = comm + ' -beam gaussASD 10 ' 




plotAngle(line1,line2)
'''
os.makedirs('saved', exist_ok=True)

plt.savefig('saved/TML.pdf', bbox_inches='tight')


# print difference
def printdiff(val,x,y):
    if ((x != 0) or (y != 0)):
        cdiff = math.fabs((x-y)/x) # calculated relative difference
        fr.write('\n\t\t'+val + ':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')
        if (cdiff > fdiff):
            #fr.write('\n\t\t'+val + ':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')
            return 1
    return 0


def printsearch(pol):
    fr.write('\n\n\tCrosSec-'+pol+': search diff')
    cext1,qext1,cabs1,qabs1 = extractCrosSeq(1,pol)
    cext2,qext2,cabs2,qabs2 = extractCrosSeq(2,pol)
    d = 0
    d += printdiff('Cext',cext1,cext2)
    d += printdiff('Qext',qext1,qext2)
    d += printdiff('Cabs',cabs1,cabs2)
    d += printdiff('Qabs',qabs1,qabs2)
    fr.write('\tCrosSec-'+pol+': done')
    return d


def compare(line1,line2):
    
    adda_run(1,line1)
    adda_run(2,line2)
    
    print('\nCompare:\n'+line1+'\n'+line2)
    fr.write('\n\nCompare:\n'+line1+'\n'+line2)
    dtotal = 0
    dtotal += printsearch('X')
    dtotal += printsearch('Y')
    if (dtotal == 0):
        print('\033[32m', '\nPassed', '\033[0m', sep='') #green stdout
    else:
        print('\033[31m', '\nNot passed (see bb_results.txt)', '\033[0m', sep='') #red
    fr.write('\n\nDone\n___')
    
    shutil.rmtree('out1')
    shutil.rmtree('out2')


fname = "bb_results.txt"
fr = open(fname, "w")
fr.write("The comparison of equivalent cases (command lines) \nfor the scattering of Bessel beams calculated in ADDA.")
fr.write("\nDifferences are shown when |diff| > " + str(fdiff) + "\n\n___")
fr.close()

fr = open(fname, "a")

cmn = ' -sym no'

print('\n\nPlane-wave limit of LE Bessel beam')
fr.write('\n\nPlane-wave limit of LE Bessel beam')
opt1 = cmn
opt2 = cmn + ' -beam besselLE 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of LM Bessel beam')
fr.write('\n\nPlane-wave limit of LM Bessel beam')
opt1 = cmn
opt2 = cmn + ' -beam besselLM 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of CS Bessel beam')
fr.write('\n\nPlane-wave limit of CS Bessel beam')
opt1 = cmn
opt2 = cmn + ' -beam besselCS 0 0'
compare(opt1,opt2)

al1 = 2
al2 = 85

print('\n\nGeneralized and LE Bessel beams')
fr.write('\n\nGeneralized and LE Bessel beams')
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' 0 0 0 1'
opt2 = cmn + ' -beam besselLE 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' 0 0 0 1'
opt2 = cmn + ' -beam besselLE 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and LM Bessel beams')
fr.write('\n\nGeneralized and LM Bessel beams')
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' 0 1 0 0'
opt2 = cmn + ' -beam besselLM 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' 0 1 0 0'
opt2 = cmn + ' -beam besselLM 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and CS Bessel beams')
fr.write('\n\nGeneralized and CS Bessel beams')
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' 0.5 0 0 0.5'
opt2 = cmn + ' -beam besselCS 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' 0.5 0 0 0.5'
opt2 = cmn + ' -beam besselCS 2 '+str(al2)
compare(opt1,opt2)

print("\n\nGeneralized and CS' Bessel beams")
fr.write("\n\nGeneralized and CS' Bessel beams")
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' 0.5 0 0 -0.5'
opt2 = cmn + ' -beam besselCSp 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' 0.5 0 0 -0.5'
opt2 = cmn + ' -beam besselCSp 2 '+str(al2)
compare(opt1,opt2)

al1 = 3 #TEL an TML types have a bigger diff for smaller angles
al2 = 85

print('\n\nGeneralized and TEL Bessel beams')
fr.write('\n\nGeneralized and TEL Bessel beams')
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' '+str(-1/math.sin(al1*math.pi/180))+' 0 0 '+str(1/math.tan(al1*math.pi/180))
opt2 = cmn + ' -beam besselTEL 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' '+str(-1/math.sin(al2*math.pi/180))+' 0 0 '+str(1/math.tan(al2*math.pi/180))
opt2 = cmn + ' -beam besselTEL 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and TML Bessel beams')
fr.write('\n\nGeneralized and TML Bessel beams')
opt1 = cmn + ' -beam besselM 2 '+str(al1)+' 0 '+str(1/math.tan(al1*math.pi/180))+' '+str(1/math.sin(al1*math.pi/180))+' 0'
opt2 = cmn + ' -beam besselTML 2 '+str(al1)
compare(opt1,opt2)

opt1 = cmn + ' -beam besselM 2 '+str(al2)+' 0 '+str(1/math.tan(al2*math.pi/180))+' '+str(1/math.sin(al2*math.pi/180))+' 0'
opt2 = cmn + ' -beam besselTML 2 '+str(al2)
compare(opt1,opt2)

fr.close()
'''