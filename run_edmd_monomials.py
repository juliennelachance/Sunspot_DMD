# Imports: 
import numpy
from math import *
import pylab
import sys
import matplotlib.colors
import scipy.linalg
import scipy.io as sio
import modred as mr
import matplotlib.pyplot as plt
from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt
from future.builtins import range
import scipy.ndimage.filters

import warnings
warnings.filterwarnings('ignore')

# Dependencies: 
# -----------------------------------------------------
# Function to sort eigenvalues and normalize for POD
def sorted_eig(M, left=False):
    w,v = scipy.linalg.eig(M, left=left, right=(not left))
    v = v.T
    sorted_wv = sorted([(abs(w[i]), list(v[i]/numpy.linalg.norm(v[i])), w[i]) for i in range(len(w))], reverse=True)
    _,v,w = zip(*sorted_wv)
    return numpy.array(w), numpy.array(v).T

# a scalar times a vector returns a vector
def scale_vector(scale, vector):
  result = [0]*len(vector)
  for i in range(len(result)):
    result[i] = scale * vector[i]
  return result

# return real part of a vector
def real_vector(vector):
  return map(lambda x: -1.0*x.real, vector)

# return imaginary part of a vector
def imag_vector(vector):
  return map(lambda x: x.imag, vector)


# Sort by closest mag to 1:
def disp_closest_evals(evals, disp_num):
	inds = numpy.argsort(numpy.abs(evals - 1.0))
	sorted_evals = evals[inds]
	print('The closest eigenvals to 1: ')
	print(sorted_evals[:disp_num])
	sorted_evals=numpy.log(sorted_evals[:disp_num])
	print('The log of these: ')
	print(sorted_evals[:disp_num])
	#print('The corresponding period in Carrington rot. number: ')
	#for x in sorted_evals[:disp_num]:  
	#	print('% 6.2f' % float(2.0*pi/numpy.imag(x)))
	print('The corresponding period in years: ')
	for x in sorted_evals[:disp_num]:  
		print('% 6.2f' % float(2.0*pi/numpy.imag(x)*27.2753/365))
	
def on_pick(event):
        ind = event.ind
	index = numpy.take(ind, 0)
	real_eig = numpy.take(numpy.real(numpy.take(v, ind)),0)
	imag_eig = numpy.take(numpy.imag(numpy.take(v, ind)),0)
        print('Eigenvalue Index: ' + repr(index))
	print('Value: ' + repr(numpy.take(numpy.take(v, ind),0)))

	imaglogeigv = numpy.imag(numpy.log(numpy.take(numpy.take(v, ind),0)))
	per = float(2.0*pi/imaglogeigv*27.2753/365)

        print('Period in years: ' + repr(per))
	print

	# Mode plot: 
	pylab.figure()
	pylab.subplot(2, 1, 1)
	dmd_m=numpy.dot(v1[:,:n_modes], d[:n_modes,:])
	pylab.plot(numpy.linspace(90,-90,num=50), numpy.real(dmd_m[:,index]), label="Real part") #, numpy.linspace(0,49,num=50)
	pylab.plot(numpy.linspace(90,-90,num=50), numpy.imag(dmd_m[:,index]), label="Imag part")
	pylab.xlim(-90,90)
	pylab.xlabel('Sine(latitude)')
	pylab.ylabel('Amplitude')
	pylab.title('Mode:')
	pylab.legend()

	# Progressing eigenvector plot: 
	pylab.subplot(2, 1, 2)
	ans_=numpy.zeros((len(A),300))
	ans_[:,0] = d[:,index]
	for i in range(1,300):
		ans_[:,i] = numpy.dot(A,ans_[:,i-1])

	pylab.imshow(numpy.real((numpy.dot(ans_[0:n_modes,:].T,v1[:,0:n_modes].T)).T), extent=[0,300,-90,90])
	pylab.colorbar()
	pylab.show(block=False)




# -----------------------------------------------------


# Load the data and store in X: 
mat_contents = sio.loadmat('sunspot_data.mat')  # MatLab file of butterfly diagram data
X = mat_contents['data_v']
# Data source: https://solarscience.msfc.nasa.gov/greenwch/bflydata.txt

# If signal falls to zero, maybe smoothing will help: 
X=scipy.ndimage.filters.gaussian_filter(X,2)


# -----------------------------------------------------
# DMD with POD modes 
# Forms an enormous matrix of the data using outer product: 
R=numpy.zeros((len(X),len(X)))
for x in X.T:
	R+=numpy.outer(x,x)
# Sorting eigenvalues and normalizing: 
w,v1 = sorted_eig(R)
# Plot to determine reasonable truncation mode values: 
#pylab.plot(w, '-x')
#pylab.title('Sorted eigenvalue plot of X')
#pylab.xlabel('Element index')
#pylab.ylabel('eigenvalue')
#pylab.show()
n_modes = 7;

# Projecting onto lower-dimensional space:
X_ = []
for i in range(n_modes):
	u_k = v1[:,i]
	X_.append(numpy.dot(X.T,u_k))
X_ = numpy.array(X_)


# -----------------------------------------------------
# EDMD, with standard basis functions: 
# (This is similar to DMD, but with a new choice of basis functions)
#Z_ = []
#for x in X_.T:
#	z = [1] #observable
#	for i,xx in enumerate(x):
#		z.append(xx)
#		for j in range(i,len(x)):
#			z.append(xx*x[j])
#			for k in range(j,len(x)):
#				z.append(xx*x[j]*x[k]);
#				for l in range(k,len(x)):
#					z.append(xx*x[j]*x[k]*x[l]);
#
#	Z_.append(z)
#Z_ = numpy.array(Z_).T

#X1 = Z_[:,:-1]
#Y1 = Z_[:,1:]
#U,s,V = numpy.linalg.svd(X1, full_matrices=False)



Z1_ = []
for x in X_.T:
	z = []
	for i,xx in enumerate(x):
		z.append(xx)		# appends x values
	Z1_.append(z)

Z2_ = []
for x in X_.T:
	z = [] #append observables
	for i,xx in enumerate(x):
		for j in range(i,len(x)):
			z.append(xx*x[j])
	Z2_.append(z)

Z3_ = []
for x in X_.T:
	z = [1] #append observables
	for i,xx in enumerate(x):
		for j in range(i,len(x)):
			for k in range(j,len(x)):
				z.append(xx*x[j]*x[k]);
				for l in range(k,len(x)):
					z.append(xx*x[j]*x[k]*x[l]);
					#for m in range(l,len(x)):
						#z.append(xx*x[j]*x[k]*x[l]*x[m]); # this one for fifth order 
						#for n in range(m,len(x)):
							#z.append(xx*x[j]*x[k]*x[l]*x[m]*x[n]);
							#for q in range(n,len(x)):
							#	z.append(xx*x[j]*x[k]*x[l]*x[m]*x[n]*x[q]);
	Z3_.append(z)

Z_=numpy.concatenate((Z1_,Z2_,Z3_),axis=1)
Z_ = numpy.array(Z_).T

print numpy.shape(Z1_)
print numpy.shape(Z2_)
print numpy.shape(Z3_)

# Now apply standard DMD to Z_ observables matrix: 
X1 = Z_[:,:-1]
Y1 = Z_[:,1:]
U,s,V = numpy.linalg.svd(X1, full_matrices=False)

modes = len(s) #can be truncated
s = s[:modes]
S = numpy.zeros((len(V), len(U)), dtype=complex)
S[:len(s), :len(s)] = numpy.diag(1.0/s)
Xpinv=numpy.dot(V.T, numpy.dot(S,U.T))

A = numpy.dot(Y1, Xpinv)

v,d = sorted_eig(A)

# Display eigenvalues closest in mag to one:
disp_num = 11
disp_closest_evals(v, disp_num)
print 


# Plot complex eigenvalues: 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(numpy.real(v), numpy.imag(v), picker=True)
ax.set_title('Eigenvalues of A')
ax.set_xlabel('Real part of eigenvalues')
ax.set_ylabel('Imag part of eigenvalues')

# Generate numbers around the complex unit circle.
N = 128
theta = scale_vector(2*pi/N, range(N))
exp_theta = map(lambda x: exp(1j * x), theta)
real_part = real_vector(exp_theta)
imag_part = imag_vector(exp_theta)
pylab.plot(real_part,imag_part)

fig.canvas.callbacks.connect('pick_event', on_pick)
pylab.show()

sys.exit(1)

# Quadratic terms plot: 

pylab.figure()
Z2_ = []
quad_len = n_modes
quad_mat = numpy.zeros(shape=(quad_len, quad_len))
mode_count = n_modes
iter_count = 0
for i in xrange(n_modes):
	for j in xrange(mode_count):
		#quad_mat[i,j]= numpy.linalg.norm(d[:,len(Z1_[0])+iter_count])	#mag of countth mode of d
		mode=numpy.dot(v1[:,:n_modes], d[:n_modes,len(Z1_[0])+iter_count])
		quad_mat[i,j]= numpy.linalg.norm(mode)	#mag of countth mode of d
		iter_count = iter_count+1
	mode_count = mode_count-1


pylab.imshow(quad_mat, interpolation='none') 
pylab.colorbar()
pylab.title('Quadratic terms: norm of DMD modes')
pylab.xlabel('POD index')
pylab.ylabel('POD index')
pylab.show()


# Log-log plotting:
v=numpy.log(v)
pylab.scatter(numpy.real(v), numpy.imag(v))
pylab.xlim([-0.1,0.1])
pylab.ylim([-0.2,0.2])
pylab.title('Log Plot of Eigenvalues of F')
pylab.xlabel('Log(Real part of eigenvalues)')
pylab.ylabel('Log(Imag part of eigenvalues)')
pylab.show()


# -----------------------------------------------------


start_val = 125		# Chosen to give interesting starting dynamics
how_many = 400		# Visualize window that extends beyond one period

pylab.imshow(numpy.real(X[:,start_val:start_val+how_many]), interpolation='none')
pylab.colorbar()
pylab.figure()

ans_=numpy.zeros((len(Z_),how_many))
ans_[:,1] = Z_[:,start_val]
for i in range(2,how_many):
	ans_[:,i] = numpy.dot(A,ans_[:,i-1])

pylab.imshow(numpy.real((numpy.dot(ans_[0:n_modes,:].T,v1[:,0:n_modes].T)).T))
pylab.colorbar()
pylab.show()
























