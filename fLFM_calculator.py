import numpy as np

#functions
def pitch(sensor_size,N):
    P = sensor_size/N
    return P

def necessary_total_magnification(FOV_target,P):
    Mt = P/FOV_target
    return Mt

def necessary_NAobj_given_res(Lambda_emission,resolution_target,Mt,delta):
    NAobj = Lambda_emission/(2*resolution_target-4*delta/Mt)
    return NAobj

def depth_of_field(Lambda_emission,NAobj,delta,Mt):
    DOF = 2*Lambda_emission/(NAobj**2) + delta/(Mt*NAobj)
    return DOF

def objective_NA_given_DOF(DOF,Lambda_emission,delta,Mt):
    NAobj = (np.sqrt(delta**2 + 8*DOF*Lambda_emission*Mt**2) - delta)/(2*DOF*Mt)
    return NAobj

def resolution_given_NA(Lambda_emission,NAobj,delta,Mt):
    res = Lambda_emission/(2*NAobj) + 2*delta/Mt
    return res

def relay_magnification(Mt,fMLA,fMO):
    Mr = fMLA/(fMO*Mt)
    return Mr
    
def second_relay_f(beam_path,Mr):
    f2 = Mr*beam_path/(2*Mr+2)
    return f2

def first_relay_f(f2,Mr):
    f1 = f2/Mr
    return f1

def field_stop(P,f2,fMLA):
    FS = P*f2/fMLA
    return FS

def resolution(Lambda_emission,NAobj,delta,Mt):
    res = Lambda_emission/(2*NAobj) + 2*delta/Mt
    return res

def depth_of_field(Lambda_emission,NAobj,delta,Mt):
    DOF = 2*Lambda_emission/(NAobj**2) + delta/(Mt*NAobj)
    return DOF

def field_of_view(P,Mt):
    FOV = P/Mt
    return FOV

# LETS GO!!
print('targets and constraints')
N=float(input('Number of elemental images (3, 5 or 7): '))
FOV_target = float(input('target field of view in mm (eg. 3): ')) # in mm, 
DOF_target = 0.2 # in mm, 200um
Lambda_emission = 0.00051 # in nm, 510nm
size_of_neuron = 0.015 # in mm, 15um
sampling_factor = 2.5 # choose in pixels per cell, just a little higher than the Nyquist rate
resolution_target = size_of_neuron/sampling_factor# in mm, 5um is alternative
print(str(N) + ' Elemental images')
print('field of view target = ' + str(round(FOV_target,2)) + ' mm')
print('resolution target = ' + str(round(resolution_target*1000,2)) + ' um')
print('depth of field target = ' + str(DOF_target) + ' mm')
print('\n')

print('decide on camera parameters')
Mpix = float(input('Camera megapixels (eg 67): '))
delta = float(input('Pixel size in um (eg 2.6): '))
delta = delta/1000 # in mm
sensor_size =  np.sqrt(Mpix*10**6)*delta # in mm
print('\nCamera: ' + str(Mpix) + ' Mpix, ' + str(delta) + ' um pixel size, ' + str(round(sensor_size,2)) + ' mm sensor size')
print('\n')

# calculations recomended MLA pitch
P = pitch(sensor_size,N)
print('maximum MLA pitch = ' + str(round(P,3)) + ' mm')

# base NA on resolution
Mt = necessary_total_magnification(FOV_target,P)
NAobj_res=necessary_NAobj_given_res(Lambda_emission,resolution_target,Mt,delta)
DOF = depth_of_field(Lambda_emission,NAobj_res,delta,Mt)
print('necessary total magnification = x' + str(round(Mt,2)))
print('necessary objective NA based on res target of ' + str(round(resolution_target*1000,2)) + ' um = ' + str(round(NAobj_res,2)))
print('depth of field = ' + str(round(DOF,2)) + ' mm')
if NAobj_res < 0:
    #error
    print('Error: negative aperture stop...')
    print('try smaller target field of view or larger pixel size')
    exit()

# NA based on depth of field
NAobj_DOF = objective_NA_given_DOF(DOF_target,Lambda_emission,delta,Mt)
res = resolution_given_NA(Lambda_emission,NAobj_DOF,delta,Mt)
print('OR\nnecessary objective NA based on DOF target of ' + str(DOF_target) + ' mm = ' + str(round(NAobj_DOF,2)))
print('Resolution = ' + str(round(res*1000,2)) + ' um')
print('\n')

# choose MLA
print('find MLA with pitch as close to max as possible and enter pitch and fMLA')
P = float(input('Choose P in mm (=<' + str(round(P,4)) +'): '))
fMLA = float(input('Choose fMLA in mm (~10xP): '))
print('MLA pitch = ' + str(round(P,2)) + ' mm')
print('MLA focal length = ' + str(fMLA) + ' mm')

# choose objective 
print('\ndecide on SLR focal length, and beam path length')
fMO = float(input('Choose fMO in mm (for SLR this might be 20 and 50): '))
AS_chosen = 2 * fMO * np.tan(np.arcsin(NAobj_res))
if AS_chosen < 0:
    #error
    print('Error: negative aperture stop...')
    exit()

# determine relay
beam_path = float(input('\nChoose beam path length in mm (~400): '))
Mr = relay_magnification(Mt,fMLA,fMO)
f2 = second_relay_f(beam_path,Mr)
f1 = first_relay_f(f2,Mr)
FS = field_stop(P,f2,fMLA)
print('relay magnification = ' + str(round(Mr,2)))
print('first relay focal length = ' + str(round(f1,2)) + ' mm')
print('second relay focal length = ' + str(round(f2,2)) + ' mm')
print('field stop = ' + str(round(FS,2)) + ' mm')

# summary
print('\n----- Final Specs ----')
print('\nComponents')
print('Camera: ' + str(Mpix) + ' Mpix, pixel size ' + str(delta) + ' mm, sensor size ' + str(round(sensor_size,2)) + ' mm')
print('MLA: pitch = ' + str(round(P,3)) + ' mm, focal length = ' + str(round(fMLA,2)) + ' mm')
print('Objective: focal length ' + str(round(fMO,2)) + ' mm, objective aperture stop = ' + str(round(AS_chosen,2)) + ' mm')
print('Relay: f1 = ' + str(round(f1,2)) + ' mm, f2 = ' + str(round(f2,2)) + ' mm' + ', field stop = ' + str(round(FS,2)) + ' mm')

# recalculate microscope specs
res = resolution(Lambda_emission,NAobj_res,delta,Mt)
DOF = depth_of_field(Lambda_emission,NAobj_res,delta,Mt)
FOV = field_of_view(P,Mt)
print('\nSpecs')
print('number of elemental images = ' + str(N))
print('Resolution = ' + str(round(res*1000,2)) + ' um')
print('depth of field = ' + str(round(DOF,2)) + ' mm')
print('field of view = ' + str(round(FOV,2)) + ' mm')
print('\n')
print('------------------------')
print('\n)')
