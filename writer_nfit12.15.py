'''This module creates par files for form factor box from par files output by 
NFIT. It accepts qzlow and qzhigh separated by a whitespace as the arguments.'''

import sys, os

def get_ChiSquare(r):
    for line in r:
        if line.startswith('# ChiSquare'):
            return float(line.split()[2])
        
def get_lambda(r):
    for line in r:
        if line.startswith('# lambda'):
            return float(line.split()[2])
            
def get_eta(r):
    for line in r:
        if line.startswith('# eta'):
            return float(line.split()[2])
            
def get_xi(r):
    for line in r:
        if line.startswith('# xi'):
            return float(line.split()[2])

def get_Kc(r):
    for line in r:
        if line.startswith('set paramArray(Kc)'):
            return float(line.split()[2])

def get_B(r):
    for line in r:
        if line.startswith('set paramArray(B)'):
            return float(line.split()[2])

def get_Lr(r):
    for line in r:
        if line.startswith('set paramArray(Lr)'):
            return float(line.split()[2])

def get_Mz(r):
    for line in r:
        if line.startswith('set paramArray(Mz)'):
            return float(line.split()[2])

def get_D(r):
    for line in r:
        if line.startswith('set paramArray(D)'):
            return float(line.split()[2])

def get_mosaic(r):
    for line in r:
        if line.startswith('set paramArray(mosaic)'):
            return float(line.split()[2])

def get_edisp(r):
    for line in r:
        if line.startswith('set paramArray(edisp)'):
            return float(line.split()[2])

def get_beamFWHM(r):
    for line in r:
        if line.startswith('set paramArray(beamFWHM)'):
            return float(line.split()[2])

def get_sDistance(r):
    for line in r:
        if line.startswith('set paramArray(sDistance)'):
            return float(line.split()[2])

def get_bc2b(r):
    for line in r:
        if line.startswith('set paramArray(bc2b)'):
            return float(line.split()[2])

def get_wavelength(r):
    for line in r:
        if line.startswith('set paramArray(wavelength)'):
            return float(line.split()[2])

def get_pixelSize(r):
    for line in r:
        if line.startswith('set paramArray(pixelSize)'):
            return float(line.split()[2])

'''Return the value of qrzero as a float'''
def get_qxzero(r):
    for line in r:
        if line.startswith('set paramArray(qxzero)'):
            return float(line.split()[2])

def get_refractiveIndex(r):
    for line in r:
        if line.startswith('set paramArray(refractiveIndex)'):
            return float(line.split()[2])
            
'''Return the value of T as a float'''
def get_T(r):
    for line in r:
        if line.startswith('set paramArray(T)'):
            return float(line.split()[2])

'''Return a list for the dataset, consisting of five floating point elements.
list[0]: sig_back
list[1]: qrlow
list[2]: qrhigh 
list[3]: qzlow 
list[4]: qzhigh'''    
def get_dataset(r):
    for line in r:
        if line.startswith('set dataset'):
            tmp = line.split()
            # for the first and last elements, remove double-quotation marks
            dataset = [float(tmp[2][1:3]), float(tmp[3]), float(tmp[4]), \
                       float(tmp[5]), float(tmp[6][0:3])]
            return dataset

'''Return the value of iter as an int'''
def get_iter(r):
    for line in r:
        if line.startswith('set iter'):
            return int(line.split()[2])

'''Return the value of aFactor as a float'''            
def get_aFactor(r):
    for line in r:
        if line.startswith('set aFactor'):
            return float(line.split()[2])

# the length of sys.argv should be 2; otherwise, complain
if len(sys.argv) != 2:
    print "Input absolute path to .par files"
    print "for example: python writer.py ~/data/fits"
    print "if par files are located in /home/biophysicists/data/fits"
else:
    output_file = open("fitdata.xls", "w")
    firstLine = 'File,Kc,B,Lr,Mz,D,mosaic,edisp,beamFWHM,sDistance,bc2b,'
    firstLine = firstLine + 'wavelength,pixelSize,qxzero,refractiveIndex,T,'
    firstLine = firstLine + 'dataset,iteration,aFactor,ChiSquare,Lambda,eta,xi\n'
    output_file.write(firstLine)
    for root, dirs, files in os.walk(sys.argv[1]):
        for filename in files:
            if filename.endswith('.par'):
                absPath = sys.argv[1] + filename
                input_file = open(absPath, "r")
                line = input_file.readline()
                # the first line indicates which version of NFIT was used
                # if par file was not created by the correct version, skip it
                if line.startswith('# NFIT version: nfit12.15'):
                    print absPath
                else:
                    continue
                ChiSquare = get_ChiSquare(input_file)
                Lambda = get_lambda(input_file)
                eta = get_eta(input_file)
                xi = get_xi(input_file)
                Kc = get_Kc(input_file)
                B = get_B(input_file)
                Lr = get_Lr(input_file)
                Mz = get_Mz(input_file)
                D = get_D(input_file)
                mosaic = get_mosaic(input_file)
                edisp = get_edisp(input_file)
                beamFWHM = get_beamFWHM(input_file)
                sDistance = get_sDistance(input_file)
                bc2b = get_bc2b(input_file)
                wavelength = get_wavelength(input_file)
                pixelSize = get_pixelSize(input_file)
                qxzero = get_qxzero(input_file)
                refractiveIndex = get_refractiveIndex(input_file)
                T = get_T(input_file)
                dataset = get_dataset(input_file)
                iteration = get_iter(input_file)
                aFactor = get_aFactor(input_file)
                input_file.close()
                
                filebasename = os.path.splitext(filename)[0]
                
                line = '%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,\
                        %s,%d,%g,%g,%g,%g,%g\n' % (filebasename, Kc, B, Lr,\
                       Mz, D, mosaic, edisp, beamFWHM, sDistance, bc2b,\
                       wavelength, pixelSize, qxzero, refractiveIndex, T,\
                       dataset, iteration, aFactor, ChiSquare, Lambda, eta, xi)
           
                output_file.write(line)
    output_file.close()

    

    
    


