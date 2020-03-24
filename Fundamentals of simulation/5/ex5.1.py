import numpy as np
import matplotlib.pyplot as plt

#Reads a square image in 8-bit/color PPM format from the given file. Note: No checks on valid format are done.
def readImage(filename):
    f = file(filename,"rb")
    
    f.readline()
    s = f.readline()
    f.readline()
    (pixel, pixel) = [t(s) for t,s in zip((int,int),s.split())]
    
    data = np.fromfile(f,dtype=np.uint8,count = pixel*pixel*3)
    img = data.reshape((pixel,pixel,3)).astype(np.double)
    
    f.close()
    
    return img, pixel
    
#Writes a square image in 8-bit/color PPM format.
def writeImage(filename, image):
    f = file(filename,"wb")
    
    pixel = image.shape[0]
    f.writelines("P6\n%d %d\n%d\n"%(pixel, pixel, 255))
    
    image = image.astype(np.uint8)
    
    image.tofile(f)
    
    f.close()
    
##############################################################################################
    #ADDED:
##############################################################################################
#Calculates sum of pixel values for one color

def sumColor(img_matrix):
    sumofvalues = np.zeros(3)	
    for colindex in np.arange(3):
        sumofvalues[colindex] = sum(sum(img_matrix[:,:,colindex]))
    print("red = {0:010.2f}, green = {1:010.2f}, blue = {2:010.2f} \n".format(sumofvalues[0],sumofvalues[1],sumofvalues[2]))
    print(sum(sumofvalues))
    return sumofvalues
##############################################################################################   
    
img, pixel = readImage("aq-original.ppm")

##############################################################################################
    #ADDED:
##############################################################################################
#check sum of pixel values for each color channel
print("original sum of pixel values \n")
sumColor(img)
#output is of the order 10**7

##############################################################################################   

#Now we set up our desired smoothing kernel. We'll use complex number for it even though it is real. 
kernel_real = np.zeros((pixel,pixel),dtype=np.complex)

hsml = 10.

##############################################################################################
    #ADDED:
##############################################################################################

knorm = 40./7./np.pi/hsml/hsml
##############################################################################################   

#now set the values of the kernel 
for i in np.arange(pixel):
    for j in np.arange(pixel):
        
        #TODO: do something sensible here to set the real part of the kernel
        #kernel_real[i, j] = ....
        
##############################################################################################
    #ADDED:
##############################################################################################
        r = np.sqrt((i - pixel/2.)**2+(j - pixel/2.)**2)
        rh = r/hsml
        if (rh>=0) and (rh < 0.5):
            kernel_real[j,i]= knorm*(1.- 6.*(rh**2) +6.*(rh**3))
        elif (rh >= 0.5) and (rh<1):
            kernel_real[j,i]= knorm*(2.*(1-rh)**3)
        else:
            kernel_real[j,i]=0
##############################################################################################   

#Let's calculate the Fourier transform of the kernel
kernel_kspace = np.fft.fft2(kernel_real)

#further space allocations for image transforms
color_real = np.zeros((pixel,pixel),dtype=np.complex)

#we now convolve each color channel with the kernel using FFTs
for colindex in np.arange(3):
    #copy input color into complex array
    color_real[:,:].real = img[:,:,colindex]
    
    
    #forward transform
    color_kspace = np.fft.fft2(color_real)
    
    #multiply with kernel in Fourier space
    #TODO: fill in code here
    
##############################################################################################
    #ADDED:
##############################################################################################
    color_kspace*=kernel_kspace
##############################################################################################   

    #backward transform
    color_real = np.fft.ifft2(color_kspace)
    
    #copy real value of complex result back into color array
    img[:,:,colindex] = color_real.real
    writeImage("aq-smoothed.ppm", img)
##############################################################################################
    #ADDED:
##############################################################################################
print("Sum of pixel values after smoothing\n")
sumColor(img)
##############################################################################################